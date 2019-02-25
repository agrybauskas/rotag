def __init_plugin__(app=None):
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('rotag', run_plugin_gui)
