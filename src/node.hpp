#ifndef _NODE_HPP_
#define _NODE_HPP_

#include <cassert>
#include <vector>
#include <ofVec2f.h>

class Node
{
public:

    Node() : density(0.0), distributions(9), u(0), node_type(FLUID)
    {

    }

    void compute_density()
    {
        density = std::accumulate(distributions.begin(),
                                  distributions.end(), 0.0);
    }

    void compute_velocity(const std::vector<ofVec2f>& c)
    {
        compute_density();


        assert (!(density == 0));

        float one_over_rho =  1.0 / density;
        u.set(0);


        for (int i = 0; i < Q; i++)
            u += (c[i]*distributions[i]);

        u *= one_over_rho;
    }

    typedef enum _node_type
    {
        FLUID = 0,
        SOLID
    } nodeType;


    float density;
    std::vector<float> distributions;
    ofVec2f u;

    nodeType node_type;

    static constexpr int D = 2;
    static constexpr int Q = 9;

protected:
private:

};

#endif //_NODE_HPP_
