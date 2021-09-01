# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 09:17:10 2021

@author: Harvinder Singh Gill
"""
from functools import wraps
#---------------------------------
# Multiple Dispatch Decorator
#---------------------------------

def get_func_kwarg_signature(kwpairs=None, criteria=None):
    func_signature_kwvaluepairs = {key:value for key, value in kwpairs.items() if key in criteria}
    print(func_signature_kwvaluepairs)
    return func_signature_kwvaluepairs

def process_kwargs(kwargs):
    func_kwnames = tuple(kw_name for kw_name in kwargs.keys())
    # or func_kwnames = fun.__code__.co_varnames
    func_kwvalues = tuple(kw_val for kw_val in kwargs.values())
    func_kwtypes = tuple(arg.__class__ for arg in func_kwvalues)
    return func_kwnames, func_kwvalues, func_kwtypes

def extract_func_kwarg_details(afunc=None, kwargs=None, dispatch_on=None, option=None):
    # A function object or just its kwargs can be passed
    #Case: function object is passed
    if afunc!=None and kwargs==None:
        kwpairs = afunc.__kwdefaults__
        func_kwnames, func_kwvalues, func_kwtypes = process_kwargs(kwpairs)
    #Case: function kwargs are passed
    if afunc==None and kwargs!=None:
        kwpairs = kwargs
        func_kwnames, func_kwvalues, func_kwtypes = process_kwargs(kwpairs)

    if dispatch_on == "kwarg_name":
        if option == "all":
            dispatch_signature = tuple(key for key in kwpairs.keys())
        if option != "all" and option !=None:
            dispatch_signature = tuple(key for key in kwpairs.keys() if key in option)

    if dispatch_on == "kwarg_value":
        if option == "all":
            dispatch_signature = tuple(val for val in kwpairs.values())
        if option != "all" and option !=None:
            dispatch_signature = tuple(val for val in kwpairs.values() if val in option)

    if dispatch_on == "kwarg_type":
        if option == "all":
            dispatch_signature = tuple(val.__class__ for val in kwpairs.values())
        if option != "all" and option !=None:
            dispatch_signature = tuple(val.__class__ for val in kwpairs.values() if val in option)

    if dispatch_on == "kwarg_name_and_type":
        if option == "all":
            dispatch_signature = {key:val.__class__  for key, val in kwpairs.items()}
        if option != "all" and option !=None:
            dispatch_signature = {key:val.__class__  for key, val in kwpairs.items() if key in option}

    if dispatch_on == "kwarg_name_and_value":
        if option == "all":
            dispatch_signature = {key:val for key, val in kwpairs.items()}
        if option != "all" and option !=None:
            dispatch_signature = {key:val for key, val in kwpairs.items() if key in option}
    return dispatch_signature, func_kwnames, func_kwvalues, func_kwtypes

decorated_method_count = 0
function_register = {} # holds registered 'function' object 
def multiple_dispatch(dispatch_on=None, option="all"):
    """
    dispatch_on = "kwarg_name" 
    dispatch_on = "kwarg_value"
    dispatch_on = "kwarg_type"
    dispatch_on = "kwarg_name_and_type" 
    dispatch_on = "kwarg_name_and_value"
    option = "all" to use all member kwargs
    option = a partial list of kwarg members
            eg: a partial list of kwarg_names:
                ["thermal_condition", "flow", "location", "unheated_case"]

    """   
    global decorated_method_count 
    decorated_method_count  = decorated_method_count  + 1               
    def register_thefunction(func):
        func_name = func.__name__
        func_kwsignature, func_kwnames, func_kwvalues, func_kwtypes = extract_func_kwarg_details(afunc=func, kwargs=None, dispatch_on=dispatch_on, option=option)
        match = function_register.get((func_name, func)) 
        if match == None:
           function_register[(func_name, func)] = (func_kwsignature, func_kwnames, func_kwvalues, func_kwtypes)
        # print("\n\n", function_register)
        @wraps(func)
        def get_kwargs_from_function_call (*args, **kwargs):
            # print("called **kwargs", kwargs)
            # print("called *args", args)
            
            called_func_kwsignature, called_kwnames, called_kwvalues, called_kwtypes = extract_func_kwarg_details(afunc=None, kwargs=kwargs, dispatch_on=dispatch_on, option=option)
            for key, vals in function_register.items():
                if (called_func_kwsignature) in vals:
                    fxn = key[1] # 2nd part of key is function object
                    return fxn(*args, **kwargs)
            if (called_func_kwsignature, called_kwnames, called_kwvalues, called_kwtypes) not in function_register.values():
                # print("called_signature:", called_func_kwsignature)
                raise TypeError("Function arguments do not match any function definition")
                # print("Check arguments")
        return get_kwargs_from_function_call      
    return register_thefunction
