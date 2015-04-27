#pragma once
#ifndef _LBM_CONTROLLER_HPP_
#define _LBM_CONTROLLER_HPP_

#include <ofVec2f.h>
#include <ofVectorMath.h>
#include "node.hpp"

class Controler
{
public:
    Controler() : x_size(0), y_size(0), size(0), w(9), c(9), nodes(0)
    {
        init_base_vectors();
        init_weights();
        init_grid();
    }

    Controler(int x, int y) : x_size(x), y_size(y), size(x*y),
        w(9), c(9), nodes(size)
    {
        init_base_vectors();
        init_weights();
        init_grid();
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

    void setup_channel()
    {
        for (int i = 0; i < x_size; i++)
        {
            int top = i;
            int bottom = i + x_size * (y_size - 1);

            nodes[top+1].node_type = Node::SOLID;
            nodes[top].node_type = Node::SOLID;
            nodes[bottom].node_type = Node::SOLID;
            nodes[bottom-1].node_type = Node::SOLID;

            nodes[top].u.set(0.0);
            nodes[bottom].u.set(0.0);

//                        for (int d = 0; d < 9; d++)
//                        {
//                            nodes[top].distributions[d] = w[i];
//                            nodes[bottom].distributions[d] = w[i];
//                        }
        }
    }

    void setup_obstacle()
    {
        int cx = 3* x_size /4;
        int cy = y_size /3;
        for (int j = cy; j < cy + (y_size/4); j++)
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
                    relax(nodes[index],0.01);
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

    //omega = 1.0/tau
    void relax(Node& n, float visc = 0.333)
    {
        Node equilibrium;

        n.compute_velocity(c);

        float omega = 1.f/(3.f * visc + 0.5f);

        float usq = n.u.lengthSquared();

        //recuring fraction
        float un1 = 1.0 - (1.5 * usq);

        for (int i = 0; i < 9; i++)
        {
            float c_dot_u = c[i].dot(n.u);

            equilibrium.distributions[i] = w[i] * n.density *
                    (un1 + c_dot_u * (3.0 + c_dot_u*4.5));

            n.distributions[i] +=
                    omega*(equilibrium.distributions[i] - n.distributions[i]);
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

    void impose_gradient()
    {
        float higher = 1. + 0.2;
        float lower  = 1. - 0.2;

        for ( int j = 2; j < y_size-2; j++)
        {
            int inlet_index = j * x_size;
            //int outlet_index = (x_size - 1) + inlet_index;
            int outlet_index = 1 + inlet_index;

            for (int d = 0; d < 9; d++)
            {
                nodes[inlet_index].distributions[d] = w[d]*higher;
                nodes[outlet_index].distributions[d] = w[d]*lower;
            }
        }


    }

    int x_size;
    int y_size;

    int size;

    //d2q9 weights
    std::vector<float> w;
    //lattice vectors
    std::vector<ofVec2f> c;

    std::vector<Node> nodes;
private:
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
