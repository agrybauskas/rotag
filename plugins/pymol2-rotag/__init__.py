def __init_plugin__( app ):
    app.menuBar.addmenuitem( 'Plugin', 'command', label='rotag',
                             command=lambda: rotag_tk_dialog( app.root ) )

def rotag_tk_dialog( parent ): # TODO: find out how to use parent object.
    import Tkinter as tk
    import ttk

    WIDTH = 500
    HEIGHT = 600

    # Initialization.
    root = tk.Tk()

    # Canvas.
    canvas = tk.Canvas( root, height=HEIGHT, width=WIDTH )
    canvas.pack()

    # Notebook.
    notebook = ttk.Notebook( root )
    notebook.place( anchor='nw' )

    # Notebook tabs.
    scan_tab          = ttk.Frame(notebook)
    library_tab       = ttk.Frame(notebook)
    about_tab         = ttk.Frame(notebook)

    notebook.add( scan_tab, text = "Scan" )
    notebook.add( library_tab, text = "Library" )
    notebook.add( about_tab, text = "About" )
