#include "ofApp.hpp"

void ofApp::setup()
{
    lbm_controler = Controler(512,32);
    lbm_controler.setup_channel();
    //lbm_controler.setup_obstacle();
    colormap.setMapFromName("jet");
}

void ofApp::update()
{
    for (int i = 0; i < 10; i++)
    {
        lbm_controler.propagate();
        lbm_controler.impose_gradient();
    }
}

void ofApp::draw()
{
    float scalex = (float)ofGetWindowWidth()/lbm_controler.x_size;
    float scaley = (float)ofGetWindowHeight()/lbm_controler.y_size;

    float scale = scaley/scalex;

    static float max = -1000;
    for (int j = 0; j < lbm_controler.y_size; j++)
    {
        for(int i = 0; i < lbm_controler.x_size; i++)
        {
            int pos = j * lbm_controler.x_size + i;

           // lbm_controler.nodes[pos].compute_velocity(lbm_controler.c);
            const ofVec2f& v = lbm_controler.nodes[pos].u;

            float mag = v.length();

            if (mag > max)
            {
                max = mag;
                std::cout << "mag: " << mag ;
                std::cout << " scaled mag = " << (12* mag * float(NCOLORS)) << std::endl;
            }

            ofVec2f p1 (i* scalex, j*scaley);
            ofVec2f p2 (p1);
            p2.x -= scalex * v.x / mag;
            p2.y -= scaley * v.y / mag;

            ofSetColor(colormap.use(12 * mag * float(NCOLORS)));
            ofDrawArrow(p1,p2, 2);
        }
    }

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
}

void ofApp::mousePressed(ofMouseEventArgs& mouse)
{
}

void ofApp::mouseReleased(ofMouseEventArgs& mouse)
{
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
