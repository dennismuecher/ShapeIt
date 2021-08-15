#include "../Include/ShapeDialogAlpha.h"
#include "../Include/ShapeFrame.h"

ShapeDialogAlpha::ShapeDialogAlpha(const TGWindow *p, const TGWindow *main, ShapeFrame *caller_obj, UInt_t w, UInt_t h, UInt_t options)
{
    fMainDialog = new TGTransientFrame(p, main, w, h, options);
    fMainDialog->Connect("CloseWindow()", "ShapeDialogAlpha", this, "CloseWindow()");
    fMainDialog->DontCallClose(); // to avoid double deletions.
    
    // use hierarchical cleaning
    fMainDialog->SetCleanup(kDeepCleanup);
        
    TGLayoutHints *fL = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);
    TGLayoutHints *fR = new TGLayoutHints( kLHintsLeft | kLHintsExpandY,
                            2, 5, 5, 2);

    
    fGDialog = new TGGroupFrame(fMainDialog, new TGString("Transformation of Oslo Literature values"),kVerticalFrame|kRaisedFrame);
    
    for (int i = 0; i < 2; i++)
       fTransform[i] = new TGCompositeFrame(fGDialog, 1, 1, kHorizontalFrame);
    
    TGLabel *t1 = new TGLabel(fTransform[0], "gSF scale B");
    fTransform[0]->AddFrame(t1, fR);
    transBDialog = new TGNumberEntry(fTransform[0], 1.0, 9,31, TGNumberFormat::kNESReal,TGNumberFormat::kNEANonNegative,TGNumberFormat::kNELLimitMinMax,0, 99999);
    
    fTransform[0]->AddFrame(transBDialog, fR);

    TGLabel *t2 = new TGLabel(fTransform[1], " slope Alpha");
    fTransform[1]->AddFrame(t2, fR);
    
    //signals
    transBDialog->Connect("ValueSet(Long_t)", "ShapeFrame", caller_obj, "TransGraph()");
    
    transAlphaDialog = new TGNumberEntry(fTransform[1], 0.0, 9,32, TGNumberFormat::kNESReal,TGNumberFormat::kNEAAnyNumber,TGNumberFormat::kNELLimitMinMax,-100, 100);

    transAlphaDialog->Connect("ValueSet(Long_t)", "ShapeFrame", caller_obj, "TransGraph()");
    fTransform[1]->AddFrame(transAlphaDialog, fR);
    TGLabel *t3 = new TGLabel(fTransform[1], " [1/MeV]");
    
    fTransform[1]->AddFrame(t3, fR);
    
    fGDialog->AddFrame(fTransform[0], fL);
    fGDialog->AddFrame(fTransform[1], fL);
    
    fMainDialog->AddFrame(fGDialog, fL);

    fMainDialog->MapSubwindows();
    fMainDialog->Resize();

    // position relative to the parent's window
    fMainDialog->CenterOnParent();

    fMainDialog->SetWindowName("Transformation");

    fMainDialog->MapWindow();
    //gClient->WaitFor(fMain);    // otherwise canvas contextmenu does not work
}

ShapeDialogAlpha::~ShapeDialogAlpha()
{
    // Delete test dialog widgets.
    
    fMainDialog->DeleteWindow();  // deletes fMain
}

double ShapeDialogAlpha::GetAlphaTransform()
{
    return transAlphaDialog->GetNumber();
}

double ShapeDialogAlpha::GetBTransform()
{
    return transBDialog->GetNumber();
}

void ShapeDialogAlpha::CloseWindow()
{
    // Called when window is closed via the window manager.
    
    delete this;
}


