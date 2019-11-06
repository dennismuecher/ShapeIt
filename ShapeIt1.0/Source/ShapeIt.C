#include "ShapeFrame.C"
#include "ShapeSetting.C"
#include "ShapeMatrix.C"
#include "ShapeGSF.C"
#include "ShapeChi2.C"



void ShapeIt() {
    static const string path = gSystem->pwd();
	ShapeFrame* test = new ShapeFrame(gClient->GetRoot(),700,500, path);
}
