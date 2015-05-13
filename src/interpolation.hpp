#pragma once
#ifndef _INTERPOLATION_HPP_
#define _INTERPOLATION_HPP_

#include <ofVec2f.h>

class Controler;

//calculate index of the field tabel at position p
//x, y calculated indexes,
//p - object / cursor position in the window coordinates
//x_size, y_sizes - limits in the field coordinates;
void calculateIndex(int&x, int& y,
                    const ofVec2f& p,
                    const int x_size, const int y_size);

/* bilinear interpolation of velocity
 * field taken form controler at position p;
 * function returns value of interpolated velocity;
 */
ofVec2f interpolate(const Controler& controler, const ofVec2f &p);



#endif //_INTERPOLATION_HPP_
