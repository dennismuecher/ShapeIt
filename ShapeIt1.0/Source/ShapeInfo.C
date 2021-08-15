class ShapeInfo {
    
private:
    TGTransientFrame  *fMain;
    TGVerticalFrame   *fFrame1;
    TGLayoutHints     *fC;
    TGPicture         *ipic;
    TGIcon            *icon;
    
public:
    ShapeInfo(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h, string path);
    virtual ~ShapeInfo();
    
    // slots
    void CloseWindow();
    void DoClose();
};


ShapeInfo::ShapeInfo(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h, string path)
{
    fMain = new TGTransientFrame(p, main, w, h ) ;
    fMain->Connect("CloseWindow()", "ShapeInfo", this, "DoClose()");
    fMain->DontCallClose(); // to avoid double deletions.
    
    // use hierarchical cleaning
    fMain->SetCleanup(kDeepCleanup);
    
    fFrame1 = new TGVerticalFrame(fMain, 500, 428, kFixedWidth | kFixedHeight);
    fC  = new TGLayoutHints(kLHintsCenterX, 1, 1, 1, 1);
    string picPath = path + "/info.jpg";
    ipic=(TGPicture *)gClient->GetPicture(picPath.c_str());
    icon = new TGIcon(fFrame1,ipic,500,428);
    
    fFrame1->AddFrame(icon, fC);
    fMain->AddFrame(fFrame1, fC);
    fMain->MapSubwindows();
    
    TGDimension size = fMain->GetDefaultSize();
    fMain->Resize(size);
    
    fMain->SetWMSize(size.fWidth, size.fHeight);
    fMain->SetWMSizeHints(size.fWidth, size.fHeight, size.fWidth, size.fHeight, 0, 0);
    fMain->SetMWMHints(kMWMDecorAll | kMWMDecorResizeH  | kMWMDecorMaximize |
                       kMWMDecorMinimize | kMWMDecorMenu,
                       kMWMFuncAll |  kMWMFuncResize    | kMWMFuncMaximize |
                       kMWMFuncMinimize,
                       kMWMInputModeless);

    //fMain->Resize();
    
    // position relative to the parent's window
    fMain->CenterOnParent();
    
    fMain->SetWindowName("About");
    
    fMain->MapWindow();
}

ShapeInfo::~ShapeInfo()
{
    // Delete dialog.
    
    fMain->DeleteWindow();   // deletes fMain
}

void ShapeInfo::CloseWindow()
{
    // Called when window is closed via the window manager.
    
    delete this;
}

void ShapeInfo::DoClose()
{
    CloseWindow();
}
