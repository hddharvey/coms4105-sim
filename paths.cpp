#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <iostream>
#include <string>
#include <cassert>
#include <limits>
#include <complex>

using namespace std;
using namespace std::complex_literals;

// Gamma correction factor for displaying the gain map
// This isn't the same as 'gamma' (the reflection coefficient)
const double DEFAULT_SQUASH_FACTOR = 1.3;

// The default wavelength that we'll work with
const double DEFAULT_WAVELENGTH = 80; // 80 pixels is reasonable

// The default reflection coefficient (can include sign)
const double DEFAULT_REFLECTIVITY = -0.7;

// Size of the dots that we draw on top of the transmitter and receiver points
const double DOT_SIZE = 10.0;

// Screen dimension constants
const int SCREEN_WIDTH = 800;
const int SCREEN_HEIGHT = 600;

enum class Grabbed
{
    NOTHING,        // nothing is being grabbed
    TRANSMITTER,    // transmitter is being grabbed
    RECEIVER,       // receiver is being grabbed
};

enum class MouseState
{
    NOTHING,        // normal
    PRESSED_GRAB,   // just grabbed it
    RELEASED_GRAB,  // released button after grabbing
    PRESSED_DROP,   // about to drop it when next released
};

struct Vec
{
    double x, y;
    Vec() : x(0.0), y(0.0) { }
    Vec(double x, double y) : x(x), y(y) { }

    Vec operator*(double a) const { return Vec(x * a, y * a); }
    Vec operator/(double a) const { return Vec(x / a, y / a); }
    Vec operator+(Vec a) const { return Vec(x + a.x, y + a.y); }
    Vec operator-(Vec a) const { return Vec(x - a.x, y - a.y); }
    double norm() const { return sqrt(x * x + y * y); }
    double dot(Vec a) const { return x * a.x + y * a.y; } 
    Vec unit() const { return *this / norm(); } // normalised version of us
    Vec perp() const { return Vec(-y, x); } // a vector perpendicular to us
    bool operator==(Vec a) const { return x == a.x && y == a.y; }
    bool operator!=(Vec a) const { return x != a.x || y != a.y; }
};

struct System
{
    Vec tx;
    Vec rx;
    Vec bounceH;      // horizontal wall only
    Vec bounceV;      // vertical wall only
    Vec bounce2[2];   // hits both walls

    System(Vec tx, Vec rx) : tx(tx), rx(rx) 
    {
        calculate_bounces();
    }

    void calculate_bounces(); // update bounce positions based on tx and rx
};

struct Colour
{
    using Compact = uint32_t;

    float r, g, b;
    Colour() 
        : r(0.0f), g(0.0f), b(0.0f) { }
    Colour(float r, float g, float b) 
        : r(r), g(g), b(b) { }
    Colour(int r, int g, int b) 
        : r(r / 255.0f), g(g / 255.0f), b(b / 255.0f) { }

    Colour(Compact rgba)
    {
        r = ((rgba & 0xFF000000) >> 24) / 255.0f;
        g = ((rgba & 0x00FF0000) >> 16) / 255.0f;
        b = ((rgba & 0x0000FF00) >> 8) / 255.0f;
    }

    Compact compact() const
    {
        return (uint8_t(r * 0xFF) << 24) + (uint8_t(g * 0xFF) << 16)
            + (uint8_t(b * 0xFF) << 8) + uint8_t(0xFF);
    }

    static const Colour PINK;
    static const Colour GREEN;
    static const Colour BLUE;
    static const Colour YELLOW;
    static const Colour RED;
    static const Colour BLACK;
    static const Colour WHITE;
};

class ColourMap
{
private:
    Colour::Compact** points;
    int numX, numY;

public:
    ColourMap(int numX, int numY);
    void set(int x, int y, Colour c);
    void draw() const;
    int num_x() const { return numX; }
    int num_y() const { return numY; }
    ~ColourMap();
};

class GainMap
{
private:
    ColourMap map;
    Vec tx; // transmitter position
    double mapMin; // minimum of brightness scale (in dB)
    double mapMax; // ... and the maximum
    double gamma; // reflection coefficient
    double lambda; // wavelength
    double squash; // raise colours to this power

    double calculate_gain(Vec pos);
    void update_map();

public:
    GainMap(Vec tx, int numX, int numY, double g, double l, double s) 
        : tx(tx), map(numX, numY), gamma(g), lambda(l), squash(s)
    {
        update_scale(); // will update the map too
    }

    void update_scale(); // also updates the map
    void draw() const { map.draw(); }

    void get_params(double& gamma, double& lambda, double& squash) const
    {
        gamma = this->gamma;
        lambda = this->lambda;
        squash = this->squash;
    }

    void set_params(double gamma, double lambda, double squash)
    {
        this->gamma = gamma;
        this->lambda = lambda;
        this->squash = squash;
        update_scale(); // updates the map too
    }

    void set_tx_pos(Vec tx) 
    { 
        this->tx = tx; 
        update_map(); 
    }
};

const Colour Colour::PINK(255, 175, 175);
const Colour Colour::GREEN(0, 255, 0);
const Colour Colour::BLUE(0, 150, 255); // lighter blue
const Colour Colour::YELLOW(255, 255, 0);
const Colour Colour::RED(255, 0, 0); 
const Colour Colour::BLACK(0, 0, 0);
const Colour Colour::WHITE(255, 255, 255);

void System::calculate_bounces()
{
    // Calculate the bounce positions for the single-bounce paths
    bounceH.x = (rx.x * tx.y + tx.x * rx.y) / (tx.y + rx.y);
    bounceH.y = 0.0;
    bounceV.x = 0.0;
    bounceV.y = (rx.y * tx.x + tx.y * rx.x) / (tx.x + rx.x);

    // Calculate the bounce positions for the double-bounce path
    bounce2[0].x = (tx.x * rx.y - rx.x * tx.y) / (tx.y + rx.y);
    bounce2[0].y = 0.0;
    bounce2[1].x = 0.0;
    bounce2[1].y = (tx.x * rx.y - rx.x * tx.y) / (tx.x + rx.x);

    // If the above is out of bounds, we'll try a vertical bounce
    if (bounce2[0].x < 0.0 || bounce2[1].y < 0.0)
    {
        bounce2[0].x = 0.0;
        bounce2[0].y = (tx.y * rx.x - rx.y * tx.x) / (tx.x + rx.x);
        bounce2[1].x = (tx.y * rx.x - rx.y * tx.x) / (tx.y + rx.y);
        bounce2[1].y = 0.0;
    }
}

ColourMap::ColourMap(int numX, int numY) : numX(numX), numY(numY)
{
    points = new Colour::Compact*[numY];
    for (int y = 0; y < numY; ++y)
    {
        points[y] = new Colour::Compact[numX];
        for (int x = 0; x < numY; ++x)
        {
            points[y][x] = Colour::BLACK.compact();
        }
    }
}

void ColourMap::set(int x, int y, Colour c)
{
    assert(x >= 0 && y >= 0 && x < numX && y < numY);
    points[y][x] = c.compact();
}

void ColourMap::draw() const
{
    const double dx = double(SCREEN_WIDTH) / (numX - 1); // fenceposts!!
    const double dy = double(SCREEN_HEIGHT) / (numY - 1);
    for (int y = 0; y < numY - 1; ++y)
    {
        for (int x = 0; x < numX - 1; ++x)
        {
            int pts[][2] = {
                {x, y},
                {x + 1, y},
                {x + 1, y + 1},
                {x, y + 1},
            };
            glBegin(GL_POLYGON);
            for (int i = 0; i < 4; ++i)
            {
                Colour c(points[pts[i][1]][pts[i][0]]);
                glColor3f(c.r, c.g, c.b);
                glVertex2d(pts[i][0] * dx, pts[i][1] * dy);
            }
            glEnd();
        }
    }
}

ColourMap::~ColourMap()
{
    for (int y = 0; y < numY; ++y)
    {
        delete[] points[y];
    }
    delete[] points;
}

void draw_line(Vec a, Vec b, double thickness, Colour colour)
{
    Vec perp = (a - b).perp().unit() * (thickness / 2);
    glBegin(GL_QUADS);
    glColor3f(colour.r, colour.g, colour.b);
    glVertex2d(a.x - perp.x, a.y - perp.y);
    glVertex2d(a.x + perp.x, a.y + perp.y);
    glVertex2d(b.x + perp.x, b.y + perp.y);
    glVertex2d(b.x - perp.x, b.y - perp.y);
    glEnd();
}

void draw_lines(const Vec* points, int count, double thickness, Colour colour)
{
    for (int i = 0; i < count - 1; ++i)
    {
        draw_line(points[i], points[i + 1], thickness, colour);
    }
}

void draw_circle(Vec center, double radius, Colour colour)
{
    const int numSegments = 10; // larger means circle is more accurate
    glBegin(GL_POLYGON);
    glColor3f(colour.r, colour.g, colour.b);
    for (int i = 0; i < numSegments; ++i)
    {
        double theta = 2.0 * M_PI * double(i) / double(numSegments);
        double x = radius * cos(theta) + center.x;
        double y = radius * sin(theta) + center.y;
        glVertex2d(x, y);
    }
    glEnd();
}

// Doesn't evaluate gains that are too close to the transmitter. Returns NaN
// instead in these cases (avoids blowing out the scale due to the asymptote).
double GainMap::calculate_gain(Vec rx)
{
    System sys(tx, rx);

    int numBounces[4] = {0, 1, 1, 2};        
    double lengths[4] = {
        (rx - tx).norm(), // direct path
        (sys.bounceH - tx).norm() + (rx - sys.bounceH).norm(),
        (sys.bounceV - tx).norm() + (rx - sys.bounceV).norm(),
        (sys.bounce2[0] - tx).norm() 
            + (sys.bounce2[1] - sys.bounce2[0]).norm()
            + (rx - sys.bounce2[1]).norm()
    };
    if (lengths[0] < DOT_SIZE / 3.0)
    {
        return numeric_limits<double>::quiet_NaN();
    }

    complex<double> gain = 0.0;
    for (int i = 0; i < 4; ++i)
    {
        double n = numBounces[i];
        double len = lengths[i];
        double beta = 2*M_PI / lambda;
        gain += (pow(gamma, n) / pow(len, 2.0)) * exp(1i*beta*len);
    }
    return 10*log10(abs(gain));
}

void GainMap::update_scale()
{
    mapMin = numeric_limits<double>::infinity();
    mapMax = -numeric_limits<double>::infinity();
    double dx = double(SCREEN_WIDTH) / (map.num_x() - 1); // fenceposts!!
    double dy = double(SCREEN_HEIGHT) / (map.num_y() - 1);
    for (int y = 0; y < map.num_y(); ++y)
    {
        for (int x = 0; x < map.num_x(); ++x)
        {
            Vec pos(x * dx, y * dy);
            double gain = calculate_gain(pos);
            if (!isnan(gain))
            {
                mapMin = min(gain, mapMin);
                mapMax = max(gain, mapMax);
            }
        }
    }
    cout << "scale: min=" << mapMin << " max=" << mapMax << endl;
    update_map();
}

void GainMap::update_map()
{
    double dx = double(SCREEN_WIDTH) / (map.num_x() - 1); // fenceposts!!
    double dy = double(SCREEN_HEIGHT) / (map.num_y() - 1);
    for (int y = 0; y < map.num_y(); ++y)
    {
        for (int x = 0; x < map.num_x(); ++x)
        {
            Vec pos(x * dx, y * dy);
            double gain = calculate_gain(pos);
            Colour c; // defaults to black
            if (isnan(gain))
            {
                gain = mapMax;
            }
            if (gain > mapMax)
            {
                c = Colour::WHITE;
            }
            else
            {
                c.r = max(0.0, (gain - mapMin) / (mapMax - mapMin));
            }
            c.r = pow(c.r, squash); // basically just gamma correction
            map.set(x, y, c);
        }
    }
}

void draw(const System& s, const GainMap& map)
{
    // Fill background of the screen
    //glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    //glClear(GL_COLOR_BUFFER_BIT);
    map.draw(); // will fill background for us

    // Draw the paths for each bounce (as well as directly from T to R)
    Vec direct[] = {s.tx, s.rx};
    Vec horizBounce[] = {s.tx, s.bounceH, s.rx};
    Vec vertBounce[] = {s.tx, s.bounceV, s.rx};
    Vec doubleBounce[] = {s.tx, s.bounce2[0], s.bounce2[1], s.rx}; 

    double thickness = 5.0;
    draw_lines(direct, 2, thickness, Colour::PINK);
    draw_lines(horizBounce, 3, thickness, Colour::GREEN);
    draw_lines(vertBounce, 3, thickness, Colour::GREEN);
    draw_lines(doubleBounce, 4, thickness, Colour::BLUE);

    // Draw the dots used to mark T and R
    draw_circle(s.tx, DOT_SIZE, Colour::YELLOW);
    draw_circle(s.rx, DOT_SIZE, Colour::RED);
}

void on_mouse_down(const SDL_MouseButtonEvent* e, Vec& tx, Vec& rx, 
    Grabbed& grabbed, MouseState& mouseState)
{
    if (e->button != SDL_BUTTON_LEFT)
    {
        return;
    }
    if (mouseState == MouseState::NOTHING)
    {
        assert(grabbed == Grabbed::NOTHING);
        double x = max(0, min(e->x, SCREEN_WIDTH));
        double y = max(0, min(e->y, SCREEN_HEIGHT));
        Vec mouse(x, SCREEN_HEIGHT - y);
        double txDistance = (mouse - tx).norm();
        double rxDistance = (mouse - rx).norm();
        cerr << "dtx:" << txDistance << "    drx:" << rxDistance << endl;
        // grabbb
        const double threshold = 50.0;
        if (txDistance < threshold)
        {
            cerr << "grabbed transmitter" << endl;
            grabbed = Grabbed::TRANSMITTER;
        }
        if (rxDistance < threshold && rxDistance < txDistance)
        {
            cerr << "grabbed receiver" << endl;
            grabbed = Grabbed::RECEIVER;
        }
        if (grabbed != Grabbed::NOTHING)
        {
            tx = (grabbed == Grabbed::TRANSMITTER) ? mouse : tx;
            rx = (grabbed == Grabbed::RECEIVER) ? mouse : rx;
            mouseState = MouseState::PRESSED_GRAB;
        }
    }
    else if (mouseState == MouseState::RELEASED_GRAB)
    {
        mouseState = MouseState::PRESSED_DROP;
    }
} 

void on_mouse_up(const SDL_MouseButtonEvent* e, Grabbed& grabbed, 
    MouseState& mouseState)
{
    if (e->button != SDL_BUTTON_LEFT)
    {
        return;
    }
    if (mouseState == MouseState::PRESSED_DROP)
    {
        // If the mouse is outside of the window, then we won't let go of it
        // and we'll go back to RELEASED_GRAB so that they have to try again.
        if (e->x < 0 || e->y < 0 || e->x >= SCREEN_WIDTH 
            || e->y >= SCREEN_HEIGHT)
        {
            mouseState = MouseState::RELEASED_GRAB;
        }
        else
        {
            mouseState = MouseState::NOTHING;
            grabbed = Grabbed::NOTHING;
        }
    }
    else if (mouseState == MouseState::PRESSED_GRAB)
    {
        mouseState = MouseState::RELEASED_GRAB;
    }
}

void on_mouse_move(const SDL_MouseMotionEvent* e, Vec& tx, Vec& rx, 
    Grabbed& grabbed)
{
    if (grabbed == Grabbed::TRANSMITTER)
    {
        tx = Vec(e->x, SCREEN_HEIGHT - e->y); 
    }
    else if (grabbed == Grabbed::RECEIVER)
    {
        rx = Vec(e->x, SCREEN_HEIGHT - e->y);
    }
}

void on_key_down(const SDL_KeyboardEvent* e, GainMap& map)
{
    if (e->keysym.sym == SDLK_SPACE)
    {
        map.update_scale();
        return;
    }

    double gamma, lambda, squash;
    map.get_params(gamma, lambda, squash);

    switch (e->keysym.sym)
    {
        case SDLK_RIGHTBRACKET:
            gamma = max(-1.0, gamma - 0.1);
            break;
        case SDLK_LEFTBRACKET:
            gamma = min(1.0, gamma + 0.1);
            break;
        case SDLK_MINUS:
            lambda = max(0.0, lambda - 10);
            break;
        case SDLK_EQUALS:
            lambda = min(SCREEN_WIDTH / 2.0, lambda + 10);
            break;
        case SDLK_COMMA:
            squash = min(2.0, squash + 0.05);
            break;
        case SDLK_PERIOD:
            squash = max(0.0, squash - 0.05);
            break;
        default:
            return;
    }
    
    map.set_params(gamma, lambda, squash);
    cout << "reflectivity=" << gamma 
        << " wavelength=" << lambda
        << " squash=" << squash 
        << endl;
}

bool init(SDL_Window*& window, SDL_GLContext& glContext)
{
    if (SDL_Init(SDL_INIT_VIDEO) < 0)
    {
        cerr << "Error: SDL_Init: " << SDL_GetError() << endl;
        return false;
    }

    window = SDL_CreateWindow(
        "Path Viewer",
        SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED,
        SCREEN_WIDTH,
        SCREEN_HEIGHT,
        SDL_WINDOW_SHOWN | SDL_WINDOW_OPENGL);
    if (!window)
    {
        cerr << "Error: SDL_CreateWindow: " << SDL_GetError() << endl;
        SDL_Quit();
        return false;
    }

    glContext = SDL_GL_CreateContext(window);
    if (!glContext)
    {
        cerr << "Error: SDL_GL_CreateContext: " << SDL_GetError() << endl;
        SDL_DestroyWindow(window);
        SDL_Quit();
        return false;
    }
    if (SDL_GL_SetSwapInterval(0) < 0)
    {
        cerr << "Error: SDL_GL_SetSwapInterval: " << SDL_GetError() << endl;
        SDL_GL_DeleteContext(glContext);
        SDL_DestroyWindow(window);
        SDL_Quit();
        return false;
    }

    // Reset the projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // Reset the modelview matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    // Default coordinates are -1.0f to 1.0f (left to right, top to bottom)
    // Scale those to -w/2 to w/2 and -h/2 to h/2
    glScalef(2.0f / SCREEN_WIDTH, 2.0f / SCREEN_HEIGHT, 1.0f);
    // Now shift so that origin is in the bottom left
    glTranslatef(-SCREEN_WIDTH / 2.0f, -SCREEN_HEIGHT / 2.0f, 0.0f);
    return true;
}

void print_help()
{
    cout << "The transmitter is yellow." << endl;
    cout << "The receiver is red." << endl;
    cout << "Click to grab either and move them around" << endl;
    cout << "Press SPACE to recalculate the scale." << endl;
    cout << "Press '-' or '+' to adjust the wavelength." << endl;
    cout << "Press '[' or ']' to adjust the reflectivity." << endl; 
    cout << "Press '<' or '>' to adjust the squash factor." << endl;
}

int main( int argc, char* args[] )
{
    SDL_Window* window;
    SDL_GLContext glContext;
    if(!init(window, glContext))
    {
        return 1;
    }
    print_help();

    Vec t0(0.2 * SCREEN_WIDTH, 0.1 * SCREEN_HEIGHT);
    Vec r0(0.5 * SCREEN_WIDTH, 0.5 * SCREEN_HEIGHT);
    System sys(t0, r0);
    GainMap map(t0, 80, 60, 
        DEFAULT_REFLECTIVITY, 
        DEFAULT_WAVELENGTH,
        DEFAULT_SQUASH_FACTOR);
    Grabbed grabbed = Grabbed::NOTHING;
    MouseState mouseState = MouseState::NOTHING;
    for (;;)
    {
        Vec oldTx = sys.tx;
        SDL_Event e;
        if (!SDL_WaitEvent(&e))
        {
            cerr << "Error: SDL_WaitEvent: " << SDL_GetError() << endl;
            break;
        }
        if (e.type == SDL_QUIT)
        {
            break;
        }
        else if (e.type == SDL_MOUSEBUTTONDOWN)
        {
            on_mouse_down(&e.button, sys.tx, sys.rx, grabbed, mouseState);
        } 
        else if (e.type == SDL_MOUSEBUTTONUP)
        {
            on_mouse_up(&e.button, grabbed, mouseState);
        }
        else if (e.type == SDL_MOUSEMOTION)
        {
            on_mouse_move(&e.motion, sys.tx, sys.rx, grabbed);
        }
        else if (e.type == SDL_KEYDOWN)
        {
            if (e.key.keysym.sym == SDLK_ESCAPE)
            {
                break;
            }
            on_key_down(&e.key, map);
        }
        if (sys.tx != oldTx)
        {
            map.set_tx_pos(sys.tx);
        }
        sys.calculate_bounces();
        draw(sys, map);
        SDL_GL_SwapWindow(window);
    }

    SDL_GL_DeleteContext(glContext);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
