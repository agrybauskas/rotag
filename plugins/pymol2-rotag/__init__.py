def __init_plugin__( app ):
    app.menuBar.addmenuitem( 'Plugin', 'command', label='rotag',
                             command=lambda: rotag_tk_dialog( app.root ) )

def rotag_tk_dialog( parent ):
    import Tkinter as tk
    import ttk

    # Initialization.
    root = tk.Tk()

    # Main frame.
    frame = tk.Frame( root )
    frame.pack()

    # Notebook.
    notebook = ttk.Notebook( frame )
    notebook.pack()

    # Notebook tabs.
    scan_tab          = ttk.Frame(notebook)
    library_tab       = ttk.Frame(notebook)
    configuration_tab = ttk.Frame(notebook)
    about_tab         = ttk.Frame(notebook)

    notebook.add( scan_tab, text = "Scan" )
    notebook.add( library_tab, text = "Library" )
    notebook.add( configuration_tab, text = "Configuration" )
    notebook.add( about_tab, text = "About" )
