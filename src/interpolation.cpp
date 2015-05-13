#include <ofAppRunner.h>

#include "interpolation.hpp"
#include "lbmController.hpp"

//calculate index of the field tabel at position p
//x, y calculated indexes,
//p - object / cursor position in the window coordinates
//x_size, y_sizes - limits in the field coordinates;
void calculateIndex(int&x, int& y,
                    const ofVec2f& p,
                    const int x_size, const int y_size)
{
    int windowWidth = ofGetWindowWidth();
    int windowHeight = ofGetWindowHeight();

    int fieldWidth = x_size;
    int fieldHeight = y_size;

    //normalized position of the particle relative to window size;
    //p is in range (0, 1);
    ofVec2f relative_pos = ofVec2f(p.x/(float)windowWidth,
                                   p.y/(float)windowHeight);

    relative_pos.x = ofWrap(relative_pos.x, 0.001, 0.999);
    relative_pos.y = ofWrap(relative_pos.y, 0.001, 0.999);
    //return 0 if we are out of the window;
    if ( (relative_pos.x < 0 || relative_pos.x > 1)
         || (relative_pos.y < 0 || relative_pos.y > 1))
    {
         x = 0;
         y = 0;
    }

    int index_x = (int)(relative_pos.x * fieldWidth);
    int index_y = (int)(relative_pos.y * fieldHeight);

    // saftey :)
    x = MAX(0, MIN(index_x, fieldWidth-1));
    y = MAX(0, MIN(index_y, fieldHeight-1));
}

/* bilinear interpolation of velocity
 * field taken form controler at position p;
 * function returns value of interpolated velocity;
 */
ofVec2f interpolate(const Controler& controler, const ofVec2f &p)
{
    ofVec2f value(0);

    int windowWidth = ofGetWindowWidth();
    int windowHeight = ofGetWindowHeight();

    int fieldWidth = controler.x_size;
    int fieldHeight = controler.y_size;

    int index_x = 0;
    int index_y = 0;

    calculateIndex(index_x, index_y, p, controler.x_size, controler.y_size);

//    //wrapping;
    int index_x2 = (index_x+1) % fieldWidth;
    int index_y2 = (index_y+1) % fieldHeight;

    index_x2 = MAX(0, MIN(index_x2, fieldWidth-1));
    index_y2 = MAX(0, MIN(index_y2, fieldHeight-1));

    ofVec2f tileSize((float)windowWidth/fieldWidth,
                     (float)windowHeight/fieldHeight);

    //velocity vertex pos;
    ofVec2f p11 ( index_x  * tileSize.x, index_y  * tileSize.y);
    ofVec2f p21 ( index_x2 * tileSize.x, index_y  * tileSize.y);
    ofVec2f p12 ( index_x  * tileSize.x, index_y2 * tileSize.y);
    ofVec2f p22 ( index_x2 * tileSize.x, index_y2 * tileSize.y);



    // indexes of the field array corresponding to verices;
    int index11 = index_y  * fieldWidth + index_x;
    int index21 = index_y  * fieldWidth + index_x2;
    int index12 = index_y2 * fieldWidth + index_x;
    int index22 = index_y2 * fieldWidth + index_x2;

    //values of the fields;
    ofVec2f value11 = controler.nodes[index11].u;
    ofVec2f value21 = controler.nodes[index21].u;
    ofVec2f value12 = controler.nodes[index12].u;
    ofVec2f value22 = controler.nodes[index22].u;

//    ofVec2f value11 = data[index11];
//    ofVec2f value21 = data[index21];
//    ofVec2f value12 = data[index12];
//    ofVec2f value22 = data[index22];

    //velocity x coordinates;
    float fxy1 = (value11.x * (p21.x - p.x) / (p21.x - p11.x)) +
                    (value21.x * (p.x - p11.x) / (p21.x - p11.x));

    float fxy2 = (value12.x * (p21.x - p.x)/(p21.x - p11.x)) +
                    (value22.x * (p.x - p11.x) / (p21.x - p11.x));

    float vx =  (fxy1 * (p12.y - p.y) / (p12.y - p11.y)) +
                    (fxy2 * (p.y - p11.y)/(p12.y - p21.y));


    //velocity y coordinate;
    fxy1 = (value11.y * (p21.x - p.x) / (p21.x - p11.x)) +
            (value21.y * (p.x - p11.x) / (p21.x - p11.x));

    fxy2 = (value12.y * (p21.x - p.x)/(p21.x - p11.x)) +
            (value22.y * (p.x - p11.x) / (p21.x - p11.x));

    float vy = (fxy1 * (p12.y - p.y) / (p12.y - p11.y)) +
            (fxy2 * (p.y - p11.y)/(p12.y - p21.y));


    value.set(vx, vy);
   // std::cout << "value: " << value << std::flush;
    return value;


}
