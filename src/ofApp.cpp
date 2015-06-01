#include "ofApp.hpp"
#include "interpolation.hpp"

void ofApp::setup()
{

    ofSetBackgroundColor(ofColor::wheat);
    viscosity.addListener(this, &ofApp::viscosityChanged);
    acceleration.addListener(this, &ofApp::accelerationChanged);
    model.addListener(this, &ofApp::modelChanged);
    //reynolds.addListener(this, &ofApp::reynoldsChanged);

    gui.setup("LBM control panel");
    //min BGK ~= 0.003
    //min MRT  =0.00000001
    //gui.add(viscosity.set("viscosity", 0.00925, 0.0001, 0.01));
    gui.add(viscosity.set("viscosity", 0.00625, 0.0001, 0.01));
    gui.add(acceleration.set("acceleration", ofVec2f(0.01, 0),
                             ofVec2f(0,0), ofVec2f(0.01,0.01)));
    gui.add(model.setup("MRT", true));
    reynolds.setup(std::string("Re: 0"));
    gui.add(&reynolds);

    Lx = Ly = 120;
    lbm_controler = Controler(Lx, Ly,
                              viscosity.get(), acceleration, MRT);

    //lbm_controler.setup_channel();
    lbm_controler.setup_cavity();
    //lbm_controler.setup_obstacle();
    colormap.setMapFromName("hsv");

//    for (int i = 0; i < 30; i++){
//        Particle p;
//        p.position.set(ofRandom(10,ofGetWindowWidth()-10),ofRandom(10, ofGetWindowHeight()-10));
//        p.velocity.set(0);
//        particles.push_back(p);
//    }
}

void ofApp::update()
{
    for (int i = 0; i < 50; i++)
    {
        lbm_controler.propagate();
    }

//    for(auto& p : particles)
//    {
//        p.update(lbm_controler);
//    }

    reynoldsChanged();

}

void ofApp::viscosityChanged(float &v)
{
    this->viscosity.set(v);
    lbm_controler.setViscosity(this->viscosity.get());
}

void ofApp::accelerationChanged(ofVec2f& a)
{
    this->acceleration.set(a);
    lbm_controler.setAcceleration(this->acceleration.get());
}

void ofApp::modelChanged(bool& b)
{
   model = b;
   if (model)
   {
       model.setName("MRT");
       lbm_controler.setModel(MRT);
   }
   else
   {
       model.setName("BGK");
       lbm_controler.setModel(BGK);

   }

   lbm_controler.u_max = -1;
}

void ofApp::reynoldsChanged()
{
    float u_max = lbm_controler.get_umax();
    float visc = lbm_controler.viscosity;
    float Re = Lx * u_max / visc;

    reynolds = "Re: " + std::to_string( int(Re));
}

void ofApp::draw()
{
    lbm_controler.draw(ofGetWindowWidth(), ofGetWindowHeight(), colormap);
//    for(auto& p : particles)
//    {
//        p.draw(colormap);
//    }
    gui.draw();
}

void ofApp::exit()
{
}

void ofApp::keyPressed(ofKeyEventArgs& key)
{
}

void ofApp::keyReleased(ofKeyEventArgs& key)
{
}

void ofApp::mouseMoved(ofMouseEventArgs& mouse)
{
}

void ofApp::mouseDragged(ofMouseEventArgs& mouse)
{
    int x, y = 0;
    calculateIndex(x, y, mouse, lbm_controler.x_size, lbm_controler.y_size);

    int index = x + lbm_controler.x_size * y;

    auto& node = lbm_controler.nodes[index];

    if (mouse.button == 1)
    {
       node.node_type = Node::SOLID;
       node.u.set(0);
    }

    if (mouse.button == 2)
    {
       node.node_type = Node::FLUID;
    }


}

void ofApp::mousePressed(ofMouseEventArgs& mouse)
{
}

void ofApp::mouseReleased(ofMouseEventArgs& mouse)
{
    if (mouse.button == 0)
    {
        particles.erase(particles.begin());
        Particle p;
        p.position.set(mouse.x,mouse.y);
        particles.push_back(p);

    }
}

void ofApp::windowResized(ofResizeEventArgs& window)
{
}

void ofApp::gotMessage(ofMessage message)
{
}

void ofApp::dragEvent(ofDragInfo dragged)
{
}
