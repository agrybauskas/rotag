def __init_plugin__( app ):
    app.menuBar.addmenuitem( 'Plugin', 'command', label='rotag',
                             command=lambda: rotag_tk_dialog( app.root ) )

def rotag_tk_dialog( parent ):
    try:
        import tkMessageBox  # Python 2
    except ImportError:
        import tkinter.messagebox as tkMessageBox  # Python 3

    tkMessageBox.showinfo( parent=parent, title='rotag',
                           message='Hello World' )
