import wx
import vtk
from vtk.wx.wxVTKRenderWindowInteractor import wxVTKRenderWindowInteractor
 
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
		
        self.filename=""
        self.reader = vtk.vtkSTLReader()
            
    def renderthis(self):
            # open a window and create a renderer
            ren = vtk.vtkRenderer()
            self.widget.GetRenderWindow().AddRenderer(ren)
   
            # o take the polygonal data from the vtkConeSource and
            # create a rendering for the renderer.
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(self.reader.GetOutputPort())

			#SetInputConnection(reader.GetOutputPort())
 
            # create an actor for our scene
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            # Add actor
            ren.AddActor(actor)
 
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
            
class TextPanel(wx.Panel):

    def __init__(self,parent):
    
        wx.Panel.__init__(self, parent, wx.ID_ANY)
        
        calc_text = wx.TextCtrl(self, wx.ID_ANY, style=wx.TE_RIGHT)
        layout = wx.BoxSizer(wx.HORIZONTAL)
        layout.Add(calc_text, 1)
        self.SetSizer(layout)
             
class MainFrame(wx.Frame):
    def __init__(self,parent,title):
        wx.Frame.__init__(self,parent,title=title,size=(650,600), style=wx.MINIMIZE_BOX|wx.SYSTEM_MENU|
                  wx.CAPTION|wx.CLOSE_BOX|wx.CLIP_CHILDREN|wx.RESIZE_BORDER)
        self.sp = wx.SplitterWindow(self)
        self.p1 = vtkPanel(self.sp)
        self.p2 = wx.Panel(self.sp,style=wx.SUNKEN_BORDER)
         
        self.sp.SplitHorizontally(self.p1,self.p2,470)
 
        self.statusbar = self.CreateStatusBar()
        self.statusbar.SetStatusText("Click on the Plot Button")
        
        self.makeMenuBar()
         
        self.plotbut = wx.Button(self.p2,-1,"plot", size=(40,20),pos=(10,10))
        self.plotbut.Bind(wx.EVT_BUTTON,self.plot)
		
    def makeMenuBar(self):
        """
        A menu bar is composed of menus, which are composed of menu items.
        This method builds a set of menus and binds handlers to be called
        when the menu item is selected.
        """

        # Make a file menu with Hello and Exit items
        fileMenu = wx.Menu()
        # The "\t..." syntax defines an accelerator key that also triggers
        # the same event
        openItem = fileMenu.Append(wx.ID_OPEN)
        helloItem = fileMenu.Append(-1, "&Hello...\tCtrl-H",
                "Help string shown in status bar for this menu item")
        fileMenu.AppendSeparator()
        # When using a stock ID we don't need to specify the menu item's
        # label
        exitItem = fileMenu.Append(wx.ID_EXIT)

        # Now a help menu for the about item
        helpMenu = wx.Menu()
        aboutItem = helpMenu.Append(wx.ID_ABOUT)

        # Make the menu bar and add the two menus to it. The '&' defines
        # that the next letter is the "mnemonic" for the menu item. On the
        # platforms that support it those letters are underlined and can be
        # triggered from the keyboard.
        menuBar = wx.MenuBar()
        menuBar.Append(fileMenu, "&File")
        menuBar.Append(helpMenu, "&Help")

        # Give the menu bar to the frame
        self.SetMenuBar(menuBar)

        # Finally, associate a handler function with the EVT_MENU event for
        # each of the menu items. That means that when that menu item is
        # activated then the associated handler function will be called.
        self.Bind(wx.EVT_MENU, self.OnOpen, openItem)
        self.Bind(wx.EVT_MENU, self.OnHello, helloItem)
        self.Bind(wx.EVT_MENU, self.OnExit,  exitItem)
        self.Bind(wx.EVT_MENU, self.OnAbout, aboutItem)


    def OnExit(self, event):
        """Close the frame, terminating the application."""
        self.Close(True)
        
    def OnOpen(self, event):
        openFileDialog = wx.FileDialog(self, "Open Mesh file", "", self.p1.filename,
                                       "*.stl", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
             
        if openFileDialog.ShowModal() == wx.ID_CANCEL:
            return
        filename = openFileDialog.GetPath()
        self.p1.reader.SetFileName(filename)

    def OnHello(self, event):
        """Say hello to the user."""
        wx.MessageBox("Hello again from wxPython")


    def OnAbout(self, event):
        """Display an About Dialog"""
        wx.MessageBox("This is a wxPython Hello World sample",
                      "About Hello World 2",
                      wx.OK|wx.ICON_INFORMATION)
 
         
 
    def plot(self,event):
        if not self.p1.isploted:
            self.p1.renderthis()
            self.statusbar.SetStatusText("Use your mouse to interact with the model")
 
         
app = wx.App(redirect=False)
frame = MainFrame(None,"FEComponent")
app.SetTopWindow(frame)
frame.Show()
app.MainLoop()