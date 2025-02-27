# -*- coding: utf-8 -*-
# Automatic Python module class tree generation, using Mermaid for rendering
# source: https://gist.github.com/Zulko/e0910cac1b27bcc3a1e6585eaee60121

import inspect

def class_name(cls):
    """Return a string representing the class"""
    # NOTE: can be changed to str(class) for more complete class info
    return cls.__name__

def classes_tree(module, base_module=None):
    if base_module is None:
        base_module == module.__name__
    module_classes = set()
    inheritances = []
    def inspect_class(cls):
        if class_name(cls) not in module_classes:
            if cls.__module__.startswith(base_module):
                module_classes.add(class_name(cls))
                for base in cls.__bases__:
                    inheritances.append((class_name(base), class_name(cls)))
                    inspect_class(base)
    for cls in [e for e in module.__dict__.values() if inspect.isclass(e)]:
        inspect_class(cls)
    return module_classes, inheritances

def classes_tree_to_mermaid(module_classes, inheritances):
    return "graph TD;\n" + "\n".join(
        list(module_classes) + [
            "%s --> %s" % (a, b)
            for a, b in inheritances
        ]
    )
if __name__ == "__main__":
  #import pizza.region as region
  #module_classes, inheritances = classes_tree(region, base_module='pizza')
  #print (classes_tree_to_mermaid(module_classes, inheritances))

    import pizza.raster as raster
    module_classes, inheritances = classes_tree(raster, base_module='pizza')
    print (classes_tree_to_mermaid(module_classes, inheritances))