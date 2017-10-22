import wx
import vtk
from vtk.wx.wxVTKRenderWindowInteractor import wxVTKRenderWindowInteractor
from vtk.util.colors import tomato
 
class vtkPanel(wx.Panel):
    def __init__(self,parent):
        wx.Panel.__init__(self, parent)
         
        #to interact with the scene using the mouse use an instance of vtkRenderWindowInteractor. 
        self.widget = wxVTKRenderWindowInteractor(self, -1)
        self.widget.Enable(1)
        self.widget.AddObserver("ExitEvent", lambda o,e,f=self: f.Close())
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.widget, 1, wx.EXPAND)
        self.SetSizer(self.sizer)
        self.Layout()
        self.isploted = False
            
    def renderthis(self):
            # open a window and create a renderer
            ren = vtk.vtkRenderer()
            self.widget.GetRenderWindow().AddRenderer(ren)
   
            # to generate polygonal data for a cone.
            cone = vtk.vtkConeSource()
            cone.SetResolution(25)
 
            # o take the polygonal data from the vtkConeSource and
            # create a rendering for the renderer.
            coneMapper = vtk.vtkPolyDataMapper()
            coneMapper.SetInput(cone.GetOutput())
 
            # create an actor for our scene
            coneActor = vtk.vtkActor()
            coneActor.SetMapper(coneMapper)
            # Add actor
            ren.AddActor(coneActor)
 
            axes = vtk.vtkAxesActor()
            self.marker = vtk.vtkOrientationMarkerWidget()
            self.marker.SetInteractor( self.widget._Iren )
            self.marker.SetOrientationMarker( axes )
            self.marker.SetViewport(0.75,0,1,0.25)
            self.marker.SetEnabled(1)
 
            ren.ResetCamera()
            ren.ResetCameraClippingRange()
            cam = ren.GetActiveCamera()
            cam.Elevation(10)
            cam.Azimuth(70)
            self.isploted = True
            
class MainMenu(wx.MenuBar):

    def __init__(self):

        wx.MenuBar.__init__(self)
        
        menu_file = wx.Menu()
        menu_file.Append(wx.ID_ANY, u"Save")
        menu_file.Append(wx.ID_ANY, u"Exit")
        menu_edit = wx.Menu()
        menu_edit.Append(wx.ID_ANY, u"Copy")
        menu_edit.Append(wx.ID_ANY, u"Paste")
        
        self.Append(menu_file, u"File")
        self.Append(menu_edit, u"Edit")
 
class MainFrame(wx.Frame):
    def __init__(self,parent,title):
        wx.Frame.__init__(self,parent,title=title,size=(650,600), style=wx.MINIMIZE_BOX|wx.SYSTEM_MENU|
                  wx.CAPTION|wx.CLOSE_BOX|wx.CLIP_CHILDREN|wx.RESIZE_BORDER)
        self.sp = wx.SplitterWindow(self,wx.ID_ANY)
        self.p1 = vtkPanel(self)
        self.p2 = wx.Panel(self.sp,style=wx.SUNKEN_BORDER)
        self.p3 = wx.TextCtrl(self.sp, -1, "This is the Edit", style=wx.TE_MULTILINE)
        
        grid_sizer = wx.FlexGridSizer(1, 2, 3, 3)
        grid_sizer.Add(self.p1,1,wx.EXPAND,0)
         
        self.sp.SplitHorizontally(self.p2,self.p3,470)
        grid_sizer.Add(self.sp,1,wx.EXPAND,0)
        self.SetSizer(grid_sizer)
        grid_sizer.Fit(self)
        grid_sizer.AddGrowableRow(0)
        grid_sizer.AddGrowableCol(0)
        grid_sizer.AddGrowableCol(1)
        self.Layout()
 
        self.statusbar = self.CreateStatusBar()
        self.statusbar.SetStatusText("Click on the Plot Button")
        
        self.SetMenuBar(MainMenu())
         
        self.plotbut = wx.Button(self.p2,-1,"plot", size=(40,20),pos=(10,10))
        self.plotbut.Bind(wx.EVT_BUTTON,self.plot)
         
 
    def plot(self,event):
        if not self.p1.isploted:
            self.p1.renderthis()
            self.statusbar.SetStatusText("Use your mouse to interact with the model")
 
         
app = wx.App(redirect=False)
frame = MainFrame(None,"FEComponent")
app.SetTopWindow(frame)
frame.Show()
app.MainLoop()