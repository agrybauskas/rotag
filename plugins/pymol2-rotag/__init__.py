def __init_plugin__(app=None):
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('rotag', run_rotag_gui)

def run_rotag_gui():
    from pymol.Qt import QtWidgets

    dialog = QtWidgets.QDialog()

    dialog.show()
