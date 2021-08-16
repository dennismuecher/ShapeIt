#ifndef SHAPEDIALOGALPHA_H
#define SHAPEDIALOGALPHA_H
//#include "ShapeFrame.h"

class ShapeFrame;
class ShapeDialogAlpha {
    
    RQ_OBJECT("ShapeDialogAlpha")
    
private:
    TGTransientFrame    *fMainDialog;
    TGGroupFrame        *fGDialog;
    TGCompositeFrame    *fTransform[2];
    TGNumberEntry       *transBDialog;              //transformation B parameter field
    TGNumberEntry       *transAlphaDialog;              //transformation alpha parameter field
   
public:
    ShapeDialogAlpha(const TGWindow *p, const TGWindow *main, ShapeFrame *caller_obj, UInt_t w, UInt_t h, double b_init = 1, double alpha_init = 0, UInt_t options = kVerticalFrame);
    
    virtual ~ShapeDialogAlpha();
    void CloseWindow();
    double GetAlphaTransform();
    double GetBTransform();   
};
#endif
