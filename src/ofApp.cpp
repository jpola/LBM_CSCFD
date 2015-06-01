#include "ofApp.hpp"
#include "interpolation.hpp"
#include <fstream>
void ofApp::setup()
{

    ofSetBackgroundColor(ofColor::wheat);
    viscosity.addListener(this, &ofApp::viscosityChanged);
    acceleration.addListener(this, &ofApp::accelerationChanged);
    model.addListener(this, &ofApp::modelChanged);
    saveBtn.addListener(this, &ofApp::saveVTK);

    //reynolds.addListener(this, &ofApp::reynoldsChanged);

    gui.setup("LBM control panel");

    //gui.add(viscosity.set("viscosity", 0.00625, 0.0001, 0.01));
    //gui.add(viscosity.set("viscosity", 0.0625, 0.0001, 0.01)); //RE 97 MRT
    //gui.add(viscosity.set("viscosity", 0.00625, 0.0001, 0.01)); //RE 970 MRT
    //gui.add(viscosity.set("viscosity", 0.00645, 0.0001, 0.01)); //RE 7747 BGK max
    gui.add(viscosity.set("viscosity", 0.0041, 0.0001, 0.01)); //RE 7747 BGK

    gui.add(acceleration.set("acceleration", ofVec2f(0.01, 0),
                             ofVec2f(0,0), ofVec2f(0.01,0.01)));
    reynolds.setup(std::string("Re: 0"));
    gui.add(&reynolds);

    gui.add(saveBtn.setup("Save to VTK"));

    gui.add(model.setup("MRT", true));
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

    for (int i = 0; i < 1550; i++)
    {
        lbm_controler.propagate();
    }

//    for(auto& p : particles)
//    {
//        p.update(lbm_controler);
//    }

    reynoldsChanged();

}

void ofApp::saveVTK()
{

    std::ofstream ofs ("simple.vtk", std::ofstream::out);
    int size = Lx * Ly;
    ofs << "# vtk DataFile Version 2.0" << std::endl;
    ofs << "LBM Cavity "  << Lx << " " << Ly << std::endl;
    ofs << "ASCII\nDATASET STRUCTURED_POINTS" << std::endl;
    ofs << "DIMENSIONS " << Lx << " " << Ly << " " << 1 << std::endl;
    ofs << "ORIGIN 0 0 0" << std::endl << "SPACING 1 1 1" << std::endl;
    ofs << "POINT_DATA " << size << std::endl;

    ofs << "VECTORS Velocity float" << std::endl;
    for (int j = Ly-1; j >= 0; j--)
    {
        for(int i = Lx-1; i >= 0; i--)
        {
            int pos = j * Lx + i;

           auto& node = lbm_controler.nodes[pos];

           if(node.node_type == Node::FLUID)
           {
               const ofVec2f& v = node.u;
               ofs << v.x << "\t"<< v.y << "\t" << 0 << std::endl;
           }
           else
           {
               ofs << 0 << "\t"<< 0 << "\t" << 0 << std::endl;
           }
        }
    }

    ofs << "VECTORS U_normalized float" << std::endl;
    float inv_unorm = 1.0 / lbm_controler.u_max;
    for (int j = Ly-1; j >= 0; j--)
    {
        for(int i = Lx-1; i >= 0; i--)
        {
            int pos = j * Lx + i;

           auto& node = lbm_controler.nodes[pos];

           if(node.node_type == Node::FLUID)
           {
               const ofVec2f& v = node.u * inv_unorm;
               ofs << v.x << "\t"<< v.y << "\t" << 0 << std::endl;
           }
           else
           {
               ofs << 0 << "\t"<< 0 << "\t" << 0 << std::endl;
           }
        }
    }
    ofs.close();
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
