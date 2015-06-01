#pragma once
#ifndef _LBM_CONTROLLER_HPP_
#define _LBM_CONTROLLER_HPP_

#include <ofVec2f.h>
#include <ofVectorMath.h>
#include "colormap/ofxColorMap.h"

#include "node.hpp"

enum LBM_MODEL {BGK, MRT};


class Controler
{
public:
    Controler() :
        x_size(0), y_size(0), size(0),
        w(9), c(9), nodes(0),
        viscosity(0), acceleration(0),
        S(9), m_eq(9), m(9), model(BGK), u_max(-100)
    {
        init_base_vectors();
        init_weights();
        init_grid();
        init_modes(viscosity, viscosity);
    }

    Controler(int x, int y,
              float viscosity, ofVec2f acceleration,
              LBM_MODEL model) :
        x_size(x), y_size(y), size(x*y),
        w(9), c(9), nodes(size),
        viscosity(viscosity), acceleration(acceleration),
        S(9), m_eq(9), m(9), model(model), u_max(-100)
    {
        init_base_vectors();
        init_weights();
        init_grid();
        init_modes(viscosity, viscosity);
    }

    void init_modes(float visc, float bulk_visc)
    {
        float delta = 0.5f;

        //pub Mrt vs. SRT
        /*
            As we mentioned before, in the LBM-MRT calculations,
            the study of Razzaghian et al. [18] is taken as a reference
            study, Thus, we will use s0 = s3 = s5 = 1 , s1 = s2 = 1.4,
            s4 = s6 = 1.2 and s7 = s8 = \omega for LBM-MRT calculations.
            And, we will compare the stability limits of the LBM-SRT
            and LBM-MRT using collision frequency ( \omega ) and 7th
            relaxation rate ( s7 = s8 ) for lid driven cavity flow,
            respectively.
         */
        S[0] = S[3] = S[5] = 1.0;
        S[2] = S[1] = 2.0 / (6.0 * bulk_visc + 1.0);
        S[7] = S[8] = 2.0 / (6.0 * visc + 1.0);

        //In case of a Poiseuille flow, Luo et al. state the following
        //formula which relates s 4 and the simulated boundary position [9]:
        //po przekształceniach
        S[4] = 1.2; //(4 - 2*S[8])/(3*delta*S[8] + 2 - S[8]);
        S[6] = S[4];
    }

    void init_grid()
    {
        for (int j = 0; j < y_size; j++)
        {
            for (int i = 0; i < x_size; i++ )
            {
                int index = i + j*x_size;
                nodes[index].node_type = Node::FLUID;
                for (int d = 0; d < 9; d++)
                {
                    nodes[index].distributions[d] = w[d];

                }
            }
        }
    }

    void setup_cavity()
    {
        for (int i = 0; i < x_size; i++)
        {
            int bottom = i + x_size * (y_size - 1);
            nodes[bottom].node_type = Node::SOLID;
            nodes[bottom].u.set(0.0);
        }

        for (int i = 0; i < y_size; i++)
        {
            //int index = x + x_size * y;
            int left  = x_size * i;
            int right = (x_size - 1) + (x_size * i);
            nodes[left].node_type = Node::SOLID;
            nodes[left].u.set(0.0);
            nodes[right].node_type = Node::SOLID;
            nodes[right].u.set(0.0);
        }

    }

    void setup_channel()
    {
        for (int i = 0; i < x_size; i++)
        {
            int top = i;
            int bottom = i + x_size * (y_size - 1);

            nodes[top].node_type = Node::SOLID;
            nodes[bottom].node_type = Node::SOLID;

            nodes[top].u.set(0.0);
            nodes[bottom].u.set(0.0);
        }
    }

    void setup_obstacle()
    {
        int cx = 3 * x_size /4;
        int cy = y_size /8;
        for (int j = 3*cy; j < 5*cy; j++)
        {
            for (int i = cx; i < cx + 1; i++ )
            {
                int index = i + x_size * j;
                nodes[index].node_type = Node::SOLID;
            }
        }
    }

    void propagate()
    {

        for (int y = 0; y < y_size; y++)
        {
            for(int x = 0; x < x_size; x++)
            {
                int index = x + y*x_size;

                if(nodes[index].node_type == Node::SOLID)
                {
                    bounce_back(nodes[index]);
                }
                else
                {
                    if (model == BGK)
                        relax(nodes[index], y);
                    else
                        relaxMrt(nodes[index], y);
                }

                int xw = (x == x_size-1    ) ? 0		: x + 1;
                int xe = (x == 0	   ) ? x_size - 1   : x - 1;
                int yn = (y == 0	   ) ? y_size - 1   : y - 1;
                int ys = (y == y_size-1    ) ? 0	    : y + 1;

                nodes[x + y*x_size].distributions[0]  = nodes[index].distributions[0];
                nodes[xe + y*x_size].distributions[1] = nodes[index].distributions[1];
                nodes[x + yn*x_size].distributions[2] = nodes[index].distributions[2];
                nodes[xw + y*x_size].distributions[3] = nodes[index].distributions[3];
                nodes[x + ys*x_size].distributions[4] = nodes[index].distributions[4];
                nodes[xe + yn*x_size].distributions[5] = nodes[index].distributions[5];
                nodes[xw + yn*x_size].distributions[6] = nodes[index].distributions[6];
                nodes[xw + ys*x_size].distributions[7] = nodes[index].distributions[7];
                nodes[xe + ys*x_size].distributions[8] = nodes[index].distributions[8];

            }
        }
    }

    void setViscosity(float viscosity)
    {
        this->viscosity = viscosity;
        init_modes(viscosity, viscosity);
    }

    void setAcceleration(ofVec2f acceleration)
    {
        this->acceleration = acceleration;
    }

    void setModel(LBM_MODEL model)
    {
        this->model = model;
    }

    void draw (int windwWidth, int windowHeight, ofxColorMap& cmap)
    {
        float scalex = (float)windwWidth/x_size;
        float scaley = (float)windowHeight/y_size;

        float scale = scaley/scalex;


        static float mag_scale = 1.0;
        for (int j = 0; j < y_size; j++)
        {
            for(int i = 0; i < x_size; i++)
            {
                int pos = j * x_size + i;

               auto& node = nodes[pos];
               ofVec2f p1 (i* scalex, j*scaley);
               ofVec2f p2 (p1);

               if(node.node_type == Node::FLUID)
               {
                   const ofVec2f& v = node.u;
                   float mag = v.length();

                   if (mag > u_max)
                   {
                       u_max = mag;
                       //std::cout << "max: " << u_max << std::endl;
                       mag_scale = float(NCOLORS) / u_max;
                   }

                   p2.x -= scalex * v.x / mag;
                   p2.y -= scaley * v.y / mag;

                   ofSetColor(cmap.use(NCOLORS - (mag * mag_scale)));
                   ofDrawArrow(p1,p2, 2);
               }
               else
               {
                   ofSetColor(ofColor::grey);
                   ofDrawRectangle(p1, scalex, scaley);

               }
            }
        }

    }

    float get_umax()
    {
        return u_max;
    }

    int x_size;
    int y_size;

    int size;

    //d2q9 weights
    std::vector<float> w;
    //lattice vectors
    std::vector<ofVec2f> c;

    std::vector<Node> nodes;

    float viscosity;

    ofVec2f acceleration;

    LBM_MODEL model;

    float u_max;


private:

    float M[9][9] = {
                        { 1, 1, 1, 1, 1, 1, 1, 1, 1},
                        {-4,-1,-1,-1,-1, 2, 2, 2, 2},
                        { 4,-2,-2,-2,-2, 1, 1, 1, 1},
                        { 0, 1, 0,-1, 0, 1,-1,-1, 1},
                        { 0,-2, 0, 2, 0, 1,-1,-1, 1},
                        { 0, 0, 1, 0,-1, 1, 1,-1,-1},
                        { 0, 0,-2, 0, 2, 1, 1,-1,-1},
                        { 0, 1,-1, 1,-1, 0, 0, 0, 0},
                        { 0, 0, 0, 0, 0, 1,-1, 1,-1}
                    };


    float M_INV [9][9]= {
        {  0.11111,  -0.11111,   0.11111,   0.00000,  -0.00000,   0.00000,  -0.00000,   0.00000,   0.00000},
        {  0.11111,  -0.02778,  -0.05556,   0.16667,  -0.16667,   0.00000,   0.00000,   0.25000,   0.00000},
        {  0.11111,  -0.02778,  -0.05556,   0.00000,   0.00000,   0.16667,  -0.16667,  -0.25000,   0.00000},
        {  0.11111,  -0.02778,  -0.05556,  -0.16667,   0.16667,   0.00000,   0.00000,   0.25000,   0.00000},
        {  0.11111,  -0.02778,  -0.05556,   0.00000,   0.00000,  -0.16667,   0.16667,  -0.25000,   0.00000},
        {  0.11111,   0.05556,   0.02778,   0.16667,   0.08333,   0.16667,   0.08333,   0.00000,   0.25000},
        {  0.11111,   0.05556,   0.02778,  -0.16667,  -0.08333,   0.16667,   0.08333,   0.00000,  -0.25000},
        {  0.11111,   0.05556,   0.02778,  -0.16667,  -0.08333,  -0.16667,  -0.08333,   0.00000,   0.25000},
        {  0.11111,   0.05556,   0.02778,   0.16667,   0.08333,  -0.16667,  -0.08333,   0.00000,  -0.25000}
    };

            /*{
                {0.1111111111111111,-0.11111111111111112,  0.1111111111111111,0,0,0,0,0,0},
                {0.1111111111111111,-0.02777777777777779, -0.055555555555555566,0.16666666666666666,-0.16666666666666669,0,0,0.25,0},
                {0.1111111111111111,-0.027777777777777762,-0.055555555555555539,0,-1.3877787807814457e-017,0.16666666666666666,-0.16666666666666669,-0.25,0},
                {0.1111111111111111,-0.02777777777777779, -0.055555555555555566,-0.16666666666666666,0.16666666666666669,0,0,0.25,0},
                {0.1111111111111111,-0.02777777777777779, -0.055555555555555566,0,1.3877787807814457e-017,-0.16666666666666666,0.16666666666666669,-0.25,0},
                {0.1111111111111111, 0.055555555555555552, 0.027777777777777776,0.16666666666666666,0.083333333333333329,0.16666666666666666,0.083333333333333329,0,0.25},
                {0.1111111111111111, 0.055555555555555552, 0.027777777777777776,-0.16666666666666666,-0.083333333333333329,0.16666666666666666,0.083333333333333329,0,-0.25},
                {0.1111111111111111, 0.055555555555555552, 0.027777777777777776,-0.16666666666666666,-0.083333333333333329,-0.16666666666666666,-0.083333333333333329,0,0.25},
                {0.1111111111111111, 0.055555555555555552, 0.027777777777777776,0.16666666666666666,0.083333333333333329,-0.16666666666666666,-0.083333333333333329,0,-0.25}
                };*/

    //MRT relaxation parameters
    // s0, s3, s5 = 0 - rho, jx, jy muszą być zachowane.
    // s4 = s6  -- related to the boundary conditions, calculated from no-slip delta = 4/3 * ( 1/s8 - 1/2)*(1/s4 - 1/2); s8 jest zadane, no-slip = delta = 1/2;
    // s1 = 2 / (6 * eta + 1)  --- eta - bulk volume viscosity, responsible for damping density and pressure fluctuations
    // s7 = s8 = 2 / (6 * nu + 1) -- nu - kinematic viscosity;
    // s2, s4, s6, reszta w granicach (0, 2)
    std::vector<float> S;
    std::vector<float> m_eq; // equilibrium moments;
    std::vector<float> m;  // node's moment value;



    //omega = 1.0/tau
    //LBM-SRT, exhibits a theoretical upper bound (omega < 2)
    void relax(Node& n, int y)
    {
        Node equilibrium;

        float tau = 3.f * viscosity + 0.5f;
        float omega = 1.f/tau;

        n.compute_velocity(c);
        if (y == 0)
            n.u += acceleration * tau;

        float usq = n.u.lengthSquared();

        for (int i = 0; i < 9; i++)
        {
            float c_dot_u = c[i].dot(n.u);

            equilibrium.distributions[i] = w[i] * n.density * ( 1. + c_dot_u * (3. + 4.5 * c_dot_u ) - 1.5 * usq);

            n.distributions[i] +=
                    omega*(equilibrium.distributions[i] - n.distributions[i]);
        }
    }

    void relaxMrt(Node& n, int y)
    {
        n.compute_velocity(c);
        if (y == 0)
        {
            float tau = 3.f * viscosity + 0.5f;
            n.u += acceleration * tau;
        }

        //calculate current value of the moments from PDF
        //mi = mi + M i,j ∗ fj
        for (int i = 0; i < 9; i++)
        {
            m[i] = 0;
            for (int j = 0; j < 9; j++)
            {
                m[i] += M[i][j] * n.distributions[j];
            }
        }

        //calculate equilibrium
        // pub MRT vs. SRT
        float rho = n.density;
        float jx = rho * n.u.x;
        float jy = rho * n.u.y;

        m_eq[0] = rho;
        m_eq[1] = (-2.0*rho) + 3.0 * (jx*jx + jy*jy); //Luo et al.
        m_eq[2] = rho - 3.0 * (jx*jx + jy*jy);
        m_eq[3] = jx;
        m_eq[4] = -jx;
        m_eq[5] = jy;
        m_eq[6] = -jy;
        m_eq[7] = jx * jx - jy * jy;
        m_eq[8] = jx*jy;

        //do collision with S as relaxation factors
        for (int i = 0; i < 9; i++)
        {
            m[i] = S[i] * (m[i] - m_eq[i]);
        }

        //update PDE by converting m -> f using M_INV;
        std::vector<float> df(9, 0);
        for (int i = 0; i < 9; i++)
        {
            df[i] = 0;
            for (int j = 0; j < 9; j++)
            {
                df[i] += M_INV[i][j] * m[j];
            }
            n.distributions[i] -= df[i];
        }

    }

    void bounce_back(Node& n)
    {
        Node temp_node = n;

        n.distributions[1] = n.distributions[3];
        n.distributions[3] = temp_node.distributions[1];

        n.distributions[2] = n.distributions[4];
        n.distributions[4] = temp_node.distributions[2];

        n.distributions[5] = n.distributions[7];
        n.distributions[7] = temp_node.distributions[5];

        n.distributions[6] = n.distributions[8];
        n.distributions[8] = temp_node.distributions[6];

    }


    void init_base_vectors()
    {
        assert(c.size() == 9);

        c[0].x =  0; c[0].y =  0;
        c[1].x =  1; c[1].y =  0;
        c[2].x =  0; c[2].y =  1;
        c[3].x = -1; c[3].y =  0;
        c[4].x =  0; c[4].y = -1;
        c[5].x =  1; c[5].y =  1;
        c[6].x = -1; c[6].y =  1;
        c[7].x = -1; c[7].y = -1;
        c[8].x =  1; c[8].y = -1;
    }

    void init_weights()
    {
        assert(w.size() == 9);

        w[0] = 4./9.;
        w[1] = w[2] = w[3] = w[4] = 1./9.;
        w[5] = w[6] = w[7] = w[8] = 1./36.;

    }
};


#endif //_LBM_CONTROLLER_HPP_
