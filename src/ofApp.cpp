#include "ofApp.hpp"
#include "interpolation.hpp"

void ofApp::setup()
{

    viscosity.addListener(this, &ofApp::viscosityChanged);
    acceleration.addListener(this, &ofApp::accelerationChanged);
    model.addListener(this, &ofApp::modelChanged);

    gui.setup("LBM control panel");
    gui.add(viscosity.set("viscosity", 0.004, 0.00001, 0.008));
    gui.add(acceleration.set("acceleration", ofVec2f(0.0005, 0),
                             ofVec2f(0,0), ofVec2f(0.002,0.002)));
    gui.add(model.setup("MRT", true));

    lbm_controler = Controler(32, 32,
                              viscosity.get(), acceleration, MRT);

    //lbm_controler.setup_channel();
    lbm_controler.setup_cavity();
    //lbm_controler.setup_obstacle();
    colormap.setMapFromName("jet");

    for (int i = 0; i < 30; i++){
        Particle p;
        p.position.set(ofRandom(10,ofGetWindowWidth()-10),ofRandom(10, ofGetWindowHeight()-10));
        p.velocity.set(0);
        particles.push_back(p);
    }
}

void ofApp::update()
{
    for (int i = 0; i < 30; i++)
    {
        lbm_controler.propagate();
    }

    for(auto& p : particles)
    {
        p.update(lbm_controler);
    }

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
}

void ofApp::draw()
{
    lbm_controler.draw(ofGetWindowWidth(), ofGetWindowHeight(), colormap);
    for(auto& p : particles)
    {
        p.draw(colormap);
    }
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
