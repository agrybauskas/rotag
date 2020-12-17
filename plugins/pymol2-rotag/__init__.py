def __init_plugin__(app):
    app.menuBar.addmenuitem('Plugin', 'command', label='rotag',
                            command=lambda: rotag_tk_dialog(app.root))

def rotag_tk_dialog(parent): # TODO: find out how to use parent object.
    import Tkinter as tk
    import ttk
    import os

    WIDTH = 500
    HEIGHT = 600

    # Initialization.
    root = tk.Tk()
    root.title("rotag")
    # root.iconbitmap(os.path.dirname(__file__) + "/icon.ico")

    # Canvas.
    canvas = tk.Canvas(root, height=HEIGHT, width=WIDTH)
    canvas.pack()

    # Notebook.
    notebook = ttk.Notebook(root)
    notebook.place(anchor='nw')

    # Notebook tabs.
    main_tab  = ttk.Frame(notebook)
    about_tab = ttk.Frame(notebook)

    # About.
    about_license = ttk.LabelFrame(about_tab, text="License")

    license_text = tk.Text(about_license)
    license_text.insert(tk.INSERT,
                        "GNU GENERAL PUBLIC LICENSE, Version 2, June 1991.",
                        tk.END)

    about_license.pack()

    notebook.add(main_tab, text="Main")
    notebook.add(about_tab, text="About")
