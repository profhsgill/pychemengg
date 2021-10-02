# -*- coding: utf-8 -*-

# This file is part of pyChemEngg python package.
 
# PyChemEngg: A python-based framework to promote problem solving and critical
# thinking in chemical engineering.

# Copyright (c) 2021 Harvinder Singh Gill <profhsgill@gmail.com>

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
r"""
Module to perform mass balances on systems without chemical reactions.

"""
from collections import Counter
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import leastsq
import logging
import copy




__all__ = ["PhysicalProcess"]

#---------------------------------
# Main outer class 
#---------------------------------

class PhysicalProcess(object):
    r"""Models a physical process such as a mixer, distillation column etc.
    
    The physical process can have streams attached to it, and a mass balance
    can be performed on these streams.
    
    
    Parameters
    ----------
    processname : `str`
        The name assigned to a physical process as its identifier.  
 
    
    Attributes
    ----------
    streamregistry : `dict` {`str`: `pychemengg.physicalmassbalance.PhysicalProcess.Stream`}
        The streamregistry is a dictionary that stores
        all streams attached to a physical process.
        The keys of the dictionary are names of streams and the 
        values are the "Stream" object itself.
        
 
    Examples
    --------
    First import the module **physicalmassbalance**
    
    .. code-block:: python
    
        >>> from pychemengg.massbalances import physicalmassbalance as pmb 
        >>> mixer1 = pmb.PhysicalProcess("mixer1")
        # This will create an instance of 'PhysicalProcess' with a name 'mixer1'.
        # Notice that here the same name is used for 'variable name' and 'process name'.
        # It does not have to be the same, however, it is less confusing if same
        # names are used.

    """  
    def __init__(self, processname):
        self.streamregistry={}
        self.processname = processname

        
    #-------------------------------------------------
    # Inner 'Stream' class and its constructor method
    #-------------------------------------------------  
    
    
    #---------------------------------------
    # 'Stream' class constructor method
    #---------------------------------------  



    def attachstreams(self, *, streamnames=None, flowrates=None, fractions=None, extrainfos=None):
        r"""Attach one or more streams to process.
        
        
        Parameters
        -----------
        streamnames : `list` [`str`, `str`, ...]
            A list containing streamnames.
        
        flowrates : `list` [`str` or `int` or `float`, `str` or `int` or `float`, ...]. optional
            A list containing flowrates for each stream. 
            
            - Use positive `int` or `float` for inlet stream.
            - Use negative `int` or `float` for outlet stream.
            
            If flowrate is unknown, use:
            
            - "F" or "+F" for inlet stream, and
            - "-F" for outlet stream
        
        fractions : `list` [ [`str` or `int` or `float`, `str` or `int` or `float`, ...], [`str` or `int` or `float`, `str` or `int` or `float`, ...], ...]. optional
            A list containing a list of component fractions for each stream.
            If the component fraction is unknown use the string "x".
        
        extrainfos : `list` [ [`str`], [`str`], ...], optional
            A list containing extrainfo for a stream represented as a string.
            
            Example:
            
            **CASE 1: Relation between component flowrates in different streams is provided.**
            
            
                Let extrainfo : ["2:S1=1.2*S4"] 
                
                This means that for component '2', the 'component flowrate' in S1
                is 1.2 times that in S4.
            
            
            **CASE 2: Relation between overall stream flowrates is provided.**
            
                Let extrainfo : ["S1=1.2*S4"]
                
                This means overall flowrate of stream S1 is 1.2 times the
                overall flowrate of stream S4.

        
        Returns
        --------
        Single `~.Stream` instance OR `list` of `~.Stream` instances
            If a single stream with just name is created, a single
            '~.Stream' instance is returned.
            If more than one stream is created, a list of `~.Stream` instances
            is returned.
        
        
        Notes
        ------
        - The Stream class requires three arguments:
            
            i) self : instance of Stream object (passed automatically by Python))
    
            ii) instance of physical process
    
            iii) streamname
    
            By calling the Stream as *PhysicalProcess.Stream(self, streamname)*,
            the arguments ii) and iii) are being provided. The 'self' in this
            statement is the physical process to which the stream is attached,
            which is argument ii) and streamname satisfies argument iii).
        
        
        - Attributes

            There are three `~.Stream` attributes set by this method.
            
            - flowrate : `int` or `float` or `str`
                Stream flowrate.
            - fractions : a list containing [`int` or `float` or `str`, `int` or `float` or `str`, ... ]
                Component mass/mole fractions of a stream.
            - extrainfo : a list of [`str`]
                Extra information between component or total flowrates between streams.
                
        
        
            +---------------------+--------------------+----------------------+
            |Method of attaching  | 'flowrate'         |'fractions' attribute |
            |'Stream'             | attribute is set as|is set as             |       
            +=====================+====================+======================+
            |Just the keyword     | "Not yet defined"  |   "Not yet defined"  |
            |"streamnames" is used|                    |                      |
            +---------------------+--------------------+----------------------+
            |The                  | User supplied      |   User supplied      |
            |keywords             | values             |   values             |
            |"streamnames" ,      |                    |                      |
            |"flowrates" ,        |                    |                      |
            |"fractions" are used |                    |                      | 
            +---------------------+--------------------+----------------------+
            
            - extrainfo attribute is either not set or set to user defined values.
            
        Examples
        ---------
        
        .. code-block:: python
        
            # First import the module **physicalmassbalance**
            >>> from pychemengg.massbalances import physicalmassbalance as pmb
            
            # Next create instance of physical process, say a mixer
            >>> mixer1 = pmb.PhysicalProcess("mixer1")
                    
            
            # How to attach a single stream with a certain name
            # -------------------------------------------------
            >>> mixer1.attachstreams(["S1"])
            # This will attach a stream named "S1" to "mixer1".
            # The returned '~.Stream' instance can be stored as follows.
            >>> S1 = mixer1.attachstreams(["S1"])
            # Flowrates and component fractions must be set
            # separately using setflow, setfractions, and setextrainfo methods.
            
            
            # How to attach many streams with given names 
            # --------------------------------------------
            >>> Water, NaOH, Product = mixer1.attachstreams(["Water", "NaOH", "Product"])
            # The returned stream objects are each stored in the stream variables
            # Water, NaOH, and Product.
            # Flowrates and component fractions must be set
            # separately using setflow, setfractions, and setextrainfo methods.
            
            
            # How to attach one or many streams with given names and data
            # -----------------------------------------------------------
            >>> streams = ["Water", "NaOH", "Product"]
            >>> flows = ["F", 1000, "-F"]
            # Use 'F' or '+F' if flowrate is unknown for inlet stream
            # Use '-F' if flowrate is unknown for outlet stream
            >>> fractions = [ ["x",0], [0,1], ["x", 0.1] ]
            # Use 'x' if a component fraction is unknown.
            >>> extrainfo = [ [], [], ["Water=0.9*Product"] ]
            >>> Water, NaOH, Product = mixer1.attachstreams(streamnames=streams,
                                                            flowrates=flows,
                                                            fractions=fractions,
                                                            extrainfos=extrainfo)
            # Here three streams are created.
            # The returned stream objects are each stored in the stream variables
            # Water, NaOH, and Product.
            # The above process can be carried for one or many streams.
        
        
        See Also
        ---------
        pychemengg.massbalances.physicalmassbalance.PhysicalProcess.Stream
        
        
        """
        stream_instances_created=[]
        
        # CASE : When a stream is created by providing its name but without data.
        
        if streamnames!=None and flowrates==None and fractions==None and extrainfos==None:
            for streamname in streamnames:
                tempstream = PhysicalProcess.Stream(self, streamname)
                tempstream.setflow("Not yet defined")
                tempstream.setfractions(["Not yet defined"])
                stream_instances_created.append(tempstream)
            if len(stream_instances_created)==1:
                # When there is a single stream, return the stream object.
                # This is done because user expects an object not a list.
                # If the list is returned, the user must do the following
                # to get the Stream object
                # somename, = mixer1.attachstreams(["somename"])
                #         | notice this comma in the line above.
                # Without the comma, the user gets a list containing the ~.Stream
                return stream_instances_created[0]
            else:
                # When there are more than one streams, return the list of stream objects.
                # In this case if the user will do the following, the user gets
                # the ~.Stream objects.
                # somename, anothername = mixer1.attachstreams(["somename", "anothername"])
                return stream_instances_created
        
        # CASE : When a stream is created by providing name and data.      
        if streamnames!=None and flowrates!=None and fractions!=None:
            # If extrainfos are not given for any stream, then
            # last argument when calling the function 'attachstreams' will be "None".
            # But a zip function is being used below, so extrainfo must be a list of same
            # length as the number of streams.
            # Therefore we transform the empty list before calling zip function
            if extrainfos==None:
                extrainfos = [[] for i in range(len(streamnames))]   
            for stream, flow, fracts, einfo in zip(streamnames,flowrates,fractions,extrainfos):
                tempstream = PhysicalProcess.Stream(self, stream)
                # If "+F" is entered convert it to "F"
                # This allows for a simpler logic during subsequent analysis.
                if flow =="+F": flow = "F"
                tempstream.setflow(flow)
                tempstream.setfractions(fracts)
                if len(einfo) != 0: # if no extrainfo for stream, 
                                    # then do not create this attribute for the stream
                    tempstream.setextrainfo(einfo)
                stream_instances_created.append(tempstream)
            
            # The stream objects are returned to caller
            # (i.e., user who is using the module). 
            # The user may capture them if they desire
            # or can simply create streams and not bother to capture
            # them in variables because the streams are anyway registered
            # in the streamregistry and are available to the code internally.
            if len(stream_instances_created)==1:
                # When there is a single stream, return the stream object.
                # This is done because user expects an object not a list.
                # If the list is returned, the user must do the following
                # to capture the Stream object, which is not intuitive
                # somename, = mixer1.attachstreams(["somename"])
                #         ^ notice this comma in the line above.
                # Without the comma, the user gets a list containing the ~.Stream            
                return stream_instances_created[0]
            else:
                # When there are more than one streams, return the list of stream objects.
                # In this case if the user will do the following, the user captures
                # the ~.Stream objects.
                # somename, anothername = mixer1.attachstreams(["somename", "anothername"])
                return stream_instances_created
            
           
    #---------------------------------------
    # Stream class ... Inner class
    #---------------------------------------  
    class Stream(object):
        r"""Creates a stream to be attached to a physical process.
        
        Parameters
        ----------
        instance of physical process : `pychemengg.physicalmassbalance.PhysicalProcess`
            Physical process to which the created stream will be attached.
        streamname : `str`
            Name of stream to be created.
            
            
        Attributes
        ----------
        streamname : `str`
            Name of stream.
        connectedprocess : `pychemengg.physicalmassbalance.PhysicalProcess`
            The instance of physical process to which the stream is attached.
        flowrate : `int` or `float` or `str`
            Stream flowrate.
        fractions : a list containing [`int` or `float` or `str`, `int` or `float` or `str`, ... ]
            Component mass/mole fractions of a stream.
        extrainfo : a list of [`str`]
            Extra information between component or total flowrates between streams.
                
        
        Returns
        -------
        None
        
                
        Notes
        -----
        Creates a stream and stores it in streamregistry (a dictionary) of the
        physical process.
        
        ~.streamregistry[streamname] = Stream object
        
        
        Examples
        --------
        .. code-block:: python
        
            To create an instance of 'Stream' use the 'attachstreams' method.
        
        See Also
        --------
        pychemengg.massbalances.physicalmassbalance.PhysicalProcess.attachstreams
        pychemengg.massbalances.physicalmassbalance.PhysicalProcess.Stream.setflow
        pychemengg.massbalances.physicalmassbalance.PhysicalProcess.Stream.setfractions
        pychemengg.massbalances.physicalmassbalance.PhysicalProcess.Stream.setextrainfo
    
    """ 
        def __init__(self, instanceof_physical_process, streamname):
            # here by Python-rules, self is instance of class 'Stream'            
            self.streamname = streamname
            self.connectedprocess = instanceof_physical_process
            self.connectedprocess.streamregistry[streamname] = self  
            # self on RHS is instance of class 'Stream'

        
        # Method of inner 'Stream' class
        #---------------------------------
        def setflow(self, flow): # setflow for a stream
            r"""Sets flowrate of a stream.
            
            
            Parameters
            ----------
            flow : `int or float or str`
                Flowrate of the stream.
                
                - Use positive `int` or `float` for inlet stream.
                - Use negative `int` or `float` for outlet stream.
                
                If flowrate is unknown, use:
                
                - "F" or "+F" for inlet stream, and
                - "-F" for outlet stream
                 
                                    
            
            Returns
            -------
            None
                
            
            
            Notes
            -----
            Creates an attribute 'flowrate' for the stream and assigns
            it a value.
                         
            
            Examples
            --------
            .. code-block:: python
            
                # For a physical process 'mixer1' and say a stream 'S1',
                # flowrate can be set as follows.
                # First import the module **physicalmassbalance**
                >>> from pychemengg.massbalances import physicalmassbalance as pmb
                
                # Next create instance of physical process, say a mixer
                >>> mixer1 = pmb.PhysicalProcess("mixer1")
                # Next, a single stream with a certain name can be attached
                # ---------------------------------------------------------
                >>> mixer1.attachstreams(["S1"])
                # This will attach a stream named "S1" to "mixer1".
                # The returned '~.Stream' instance can be stored as follows.
                >>> S1 = mixer1.attachstreams(["S1"])
                # If flowrate is 1335, it can be assigned as
                >>> S1.setflow(1335)
                >>> print(S1.flowrate)
                1335
                # If flowrate is -300 (assuming stream is outlet)
                >>> S1.setflow(-300)
                >>> print(S1.flowrate)
                -300
                
            
            
            See Also
            --------
            pychemengg.massbalances.physicalmassbalance.PhysicalProcess.attachstreams
            pychemengg.massbalances.physicalmassbalance.PhysicalProcess.Stream.setfractions
            pychemengg.massbalances.physicalmassbalance.PhysicalProcess.Stream.setextrainfo       
            
            """  
            # If "+F" is entered convert it to "F"
            # This allows for a simpler logic during subsequent analysis.
            if flow =="+F": flow = "F"
            self.flowrate = flow            

        # Method of inner 'Stream' class
        #---------------------------------
        def setfractions(self, fracts):
            r"""Set component mole/mass fractions of a stream.
            
            
            Parameters
            ----------
            fracts: A list of [`int or float or str`, `int or float or str`, ...]
                Mass or mole fractions of stream components.
                
                If a fraction for a component is unknown, use:
                
                - "x" 
                 
                                    
            
            Returns
            -------
            None
                
            
            
            Notes
            -----
            Creates an attribute 'fractions' for the stream and assigns
            it a value.
                         
            
            Examples
            --------
            .. code-block:: python
            
                # For a physical process 'mixer1' and say a stream 'S1',
                # flowrate can be set as follows.
                # First import the module **physicalmassbalance**
                >>> from pychemengg.massbalances import physicalmassbalance as pmb
                
                # Next create instance of physical process, say a mixer
                >>> mixer1 = pmb.PhysicalProcess("mixer1")
                # Next, a single stream with a certain name can be attached
                # ---------------------------------------------------------
                >>> mixer1.attachstreams(["S1"])
                # This will attach a stream named "S1" to "mixer1".
                # The returned '~.Stream' instance can be stored as follows.
                >>> S1 = mixer1.attachstreams(["S1"])
                # If the stream has four components and with mass fractions
                # 0.1, 0.3, 0.5, 0.1, then specify them as follows:
                >>> S1.setfractions([0.1, 0.3, 0.5, 0.1]
                # If one of them in an unknown, say component '3'
                >>> S.setfractions([0.1, 0.3, "x", 0.1])
                

            See Also
            --------
            pychemengg.massbalances.physicalmassbalance.PhysicalProcess.attachstreams
            pychemengg.massbalances.physicalmassbalance.PhysicalProcess.Stream.setflow
            pychemengg.massbalances.physicalmassbalance.PhysicalProcess.Stream.setextrainfo       
            
            """
            self.fractions = fracts
            
        # Method of inner 'Stream' class
        #---------------------------------
        def setextrainfo(self, info): 
            r"""Set extra information of a stream.
            
            
            Parameters
            ----------
            extrainfo: A list of [`str`, ...]
                String representation of relationship between overall or
                component flowrates amongst different streams, or
                some other relationship that can help solve the mass balance.
                                  
            
            Returns
            -------
            None
                
            
            
            Notes
            -----
            Creates an attribute 'extrainfo' for the stream and assigns
            it a value.
                         
            
            Examples
            --------
            .. code-block:: python
            
                # For a physical process 'mixer1' and say a stream 'S1',
                # flowrate can be set as follows.
                # First import the module **physicalmassbalance**
                >>> from pychemengg.massbalances import physicalmassbalance as pmb
                
                # Next create instance of physical process, say a mixer
                >>> mixer1 = pmb.PhysicalProcess("mixer1")
                # Next, a single stream with a certain name can be attached
                # ---------------------------------------------------------
                >>> mixer1.attachstreams(["S1"])
                # This will attach a stream named "S1" to "mixer1".
                # The returned '~.Stream' instance can be stored as follows.
                >>> S1 = mixer1.attachstreams(["S1"])
                # Say that stream S1 overall flowrate is half that of S4
                # then S1=0.5*S4
                >>> S1.extrainfo[["S1=0.5*S4"]]
                            
    
            See Also
            --------
            pychemengg.massbalances.physicalmassbalance.PhysicalProcess.attachstreams
            pychemengg.massbalances.physicalmassbalance.PhysicalProcess.Stream.setflow
            pychemengg.massbalances.physicalmassbalance.PhysicalProcess.Stream.setfractions      
            
            """
            self.extrainfo = info
        
        # Method of inner 'Stream' class
        #---------------------------------
        def deleteextrainfo(self): # delete extra info of a stream
            r"""Delete extra information of a stream.
            
            
            Parameters
            ----------
            None
                                  
            
            Returns
            -------
            None
                
            
            
            Notes
            -----
            Deletes 'extrainfo' attribute of a stream if it exists
                         
            
            Examples
            --------
            .. code-block:: python
            
                >>> S1.deleteextrainfo()
                            
    
            See Also
            --------
            pychemengg.massbalances.physicalmassbalance.PhysicalProcess.attachstreams
            pychemengg.massbalances.physicalmassbalance.PhysicalProcess.Stream.setflow
            pychemengg.massbalances.physicalmassbalance.PhysicalProcess.Stream.setfractions      
            
            """
            if hasattr(self, "extrainfo"):
                delattr(self, "extrainfo") 
        
        # Method of inner 'Stream' class
        #---------------------------------
        # copy one stream to another
        def setequalto(self, streamtocopy,  flowdirection="+"):
    
            r"""Copy data of one stream to another stream.
            
            
            Parameters
            ----------
            streamtocopy : `pychemengg.physicalmassbalance.PhysicalProcess.Stream`
                Stream to be copied.
            destination : `pychemengg.physicalmassbalance.PhysicalProcess.Stream`
                Destination stream. "self" is the destination stream.
            flowdirection : `str`
                "+" or "-" depending on whether destination stream is entering
                or exiting the physical process, default = "+"
                                  
            
            Returns
            -------
            None
                
            
            Examples
            --------
            .. code-block:: python
            
                # Consider two physical processes 'process1' and 'process2'.
                # Consider a stream S1 from process1 enters process2 as S2.
                # Then it would be desirable to set S2 equal to S1 so that
                # data of S1 can be copied into S2.
                # This can be done as follows.
                # First import the module **physicalmassbalance**
                >>> from pychemengg.massbalances import physicalmassbalance as pmb
                
                # Next create instance of physical process process1
                >>> process1 = pmb.PhysicalProcess("process1")

                # Attach S1 to process1.
                >>> process1.S1 = process1.attachstreams(streamnames=["S1"], flowrates=[-300], fractions=[[0.1, 0.3, 0.6]])
                
                >>> print(process1)
                # Output is as follows:
                    
                ====================================
                Process streams for 'process1' are :
                ====================================
                Stream = S1; Flowrate = -300.00; Fractions = ['0.1000', '0.3000', '0.6000']; Extra Info = None 
                ------------------------------------
                END
                
                # Next create instance of physical process process2
                >>> process2 = pmb.PhysicalProcess("process2")
                # Attach S2 to process2
                >>> process2.S2 = process2.attachstreams(streamnames=["S2"])
                >>> process2.S2.setequalto(process1.S1, flowdirection="+")
                # S2 is inlet thus it has '+' sign for flowdirection
                >>> print(process2)
                
                ====================================
                Process streams for 'process2' are :
                ====================================
                Stream = S2; Flowrate = 300.00; Fractions = ['0.1000', '0.3000', '0.6000']; Extra Info = None 
                ------------------------------------
                END

            """
            tempstream = copy.deepcopy(streamtocopy) # this is done so that we do not change original
                                                     # stream passed to the method
            if hasattr(tempstream, "flowrate"):
                directionof_stream_beingcopied = self.connectedprocess._getsignof(tempstream.flowrate)
                direction_required = eval(flowdirection + "1")
            
            # Direction of stream being coppied can be '+' or '-'
            # The following logic correctly assigns the desired flowdirection
            # to copied stream. It applied __pos__ and __neg__ methods
            # such as -tempstream or +tempstream to first get the correct sign
            # and then assign this to copied stream
            if direction_required * directionof_stream_beingcopied < 0:  
                self.setflow((-tempstream).flowrate) # here first -tempstream changes flow sign
                                                    # then this new flowrate is assigned to copied object
            if direction_required * directionof_stream_beingcopied > 0:  
                self.setflow((+tempstream).flowrate)# here first -tempstream changes flow sign
                                                    # then this new flowrate is assigned to copied object
            if hasattr(streamtocopy, "fractions"):
                self.setfractions(tempstream.fractions)
            if hasattr(streamtocopy, "extrainfo"):
                self.setextrainfo(tempstream.extrainfo)

        # Method of inner 'Stream' class
        #---------------------------------   
        def __pos__(self):
            r"""Performs unary operation '+' on flow rate of a stream.
            
            
            Parameters
            ----------
            stream : `pychemengg.physicalmassbalance.PhysicalProcess.Stream`
                Instance of `~.Stream` whose flow rate will be operated on.
                                  
            
            Returns
            -------
            stream : `pychemengg.physicalmassbalance.PhysicalProcess.Stream`
                Instance of `~.Stream` without change to flow rate is returned.
                
            
            Examples
            --------
            .. code-block:: python
            
                # For a stream 'S1' attached to 'process1'
                >>> +S1
                # The stream is returned without change to flowrate

            """    
        
            if self.flowrate == "F":
                self.setflow("F")
                return self
            if self.flowrate == "-F":
                self.setflow("-F")
                return self
            if (isinstance(self.flowrate, int) or isinstance(self.flowrate, float)):  
                self.setflow(+self.flowrate)
                return self
        

        # Method of inner 'Stream' class
        #---------------------------------   
        def __neg__(self): # unary operator '-' to negate flow rate of a passed stream
            r"""Performs unary operation '-' on flow rate of a stream.
            
            
            Parameters
            ----------
            stream : `pychemengg.physicalmassbalance.PhysicalProcess.Stream`
                Instance of `~.Stream` whose flow rate will be operated on.
                                  
            
            Returns
            -------
            stream : `pychemengg.physicalmassbalance.PhysicalProcess.Stream`
                Instance of `~.Stream` with negated flow rate is returned.
                
            
            Examples
            --------
            .. code-block:: python
            
                # For a stream 'S1' attached to 'process1'
                >>> +S1
                # The stream is returned with negated flowrate
            
            """
                
            if self.flowrate == "F":
                self.setflow("-F")
                return self
            if self.flowrate == "-F":
                self.setflow("F")
                return self
            if (isinstance(self.flowrate, int) or isinstance(self.flowrate, float)):  
                self.setflow(-self.flowrate)
                return self

        # Method of inner 'Stream' class
        #---------------------------------        
        def __add__(self, other): # overload binary '+' operator on object 'Stream'
            
            r""" Method to add two streams
            1) flowrates get added by first applying 'abs' on flowrates
            2) component fractions are computed for the mixed stream
            """
            tempprocess = PhysicalProcess("tempprocess") # create instance of class
            tempstream = tempprocess.attachstreams(streamnames=["temp"]) # create a temporary stream to return
            try: 
                totalflow = abs(self.flowrate) + abs(other.flowrate)
            except Exception:
                print("\nUnable to add stream flowrates.")
                print("Please check flow rate values to ensure they are all known.")
                logging.exception("Caught an error in method 'add' of Stream class while computing totalflow")
                return
            else:
                tempstream.setflow(totalflow) # set flow rate of temp stream
                print("\nTotal flow = ", totalflow)
            
            # compute component fractions of mixed stream
            tempfractions = []
            for x, y in zip(self.fractions, other.fractions):
                try:
                    tempfractions.append((x*abs(self.flowrate) + y*abs(other.flowrate)) / totalflow)
                except Exception:
                    print("\nCannot compute fractions in added streams because:")
                    print("Unknown fractions or flowrates may exist in one of the streams.\n")
                    logging.exception("Caught an error in method 'add' of Stream class while computing tempfractions")
            if len(tempfractions) > 0:
                print("fractions", tempfractions)        
            tempstream.setfractions(tempfractions) # set fractions of temp stream
            del tempprocess # delete "tempprocess" 
            return tempstream # return temp stream 

        # Method of inner 'Stream' class
        #---------------------------------        
        def __sub__(self, other): # overload binary '-' operator on object 'Stream'
            r""" Method to subtract two streams
            1) flowrates get subtracted by first applying 'abs' on flowrates
            2) component fractions are computed for the mixed stream
            """
            tempprocess = PhysicalProcess("tempprocess") # create instance of class to s
            tempstream = tempprocess.attachstreams(streamnames=["temp"]) # create a temporary stream to return
            try: 
                totalflow = abs(self.flowrate) - abs(other.flowrate)
            except Exception:
                print("\nUnable to add stream flowrates.")
                print("Please check flow rate values to ensure they are all known.")
                logging.exception("Caught an error in method 'subtract' of Stream class while computing totalflow")
                return
            else:
                tempstream.setflow(totalflow) # set flow rate of temp stream
                print("\nTotal flow = ", totalflow)
            tempfractions = []
            for x, y in zip(self.fractions, other.fractions):
                try:
                    tempfractions.append((x*abs(self.flowrate) - y*abs(other.flowrate)) / totalflow)
                except Exception:
                    print("\nCannot compute fractions in added streams because:")
                    print("Unknown fractions or flowrates may exist in one of the streams.\n")
                    logging.exception("Caught an error in method 'subtract' of Stream class while computing tempfractions")
            if len(tempfractions) > 0:
                print("fractions", tempfractions)        
            tempstream.setfractions(tempfractions) # set fractions of temp stream
            del tempprocess # delete "tempprocess" 
            return tempstream # return temp stream 

        # Method of inner 'Stream' class
        #---------------------------------         
        def _calc_missing_xfraction(self):
            r"""
            Computes the missing 'x' fraction in a 
            given stream by using Σx = 1

            Returns
            -------
            List containing [`bool`, `float` or `str`]
            
            - If `bool` == True, `float`= computed 'x' value
            - If `bool` == False, `str` = text message

            """
            if self.fractions.count("x") == 1.0:
                missingx = (1 - sum([val for val in self.fractions if val != "x"])) # because Σx = 1
                return [True, missingx]
            if self.fractions.count("x") > 1.0:
                text = "More than one 'x' unknown fractions in the stream. Σx = 1 is insufficient to find them."
                return [False, text]
            if self.fractions.count("x") == 0.0:
                text = "There are no missing fractions."
                return [False, text]   

        # Method of inner 'Stream' class
        #---------------------------------         
        def __str__(self):
            r""" Implements 'print' function on Stream instance
                       
            """

            if hasattr(self,"extrainfo"):
                storedextrainfo, = self.extrainfo
            else:
                storedextrainfo = "None"
            unit_name = (self.connectedprocess.processname).upper()
            
            ##### Older version start
            # return('Stream: {streamname}\n'
            # 'Attached to the physical process: {connectedprocess.processname}\n'
            # 'flowrate: {flowrate}\n'
            # 'fractions: {fractions}\n'
            # ).format(**self.__dict__) + "extrainfo: " + storedextrainfo
            ##### Older version ends
            
            text_to_print = (f"\nStream: {self.streamname}"
                             f"\nAttached to the physical process: {unit_name}"
                             f"\nflowrate: {self.flowrate}"
                             f"\nfractions: {self.fractions}"
                             f"\nextrainfo: {storedextrainfo}")
            return text_to_print
       
    #-----------------------------------------
    # Methods for outer class "PhysicalProcess"
    #-----------------------------------------
    
    # Method for outer class "PhysicalProcess"
    #----------------------------------------    
    def _make_numpyarrayof_fractions(self):
        r""" Creates numpy array of all stream-fractions for a physical process.
                   
        """
        a = [streamobject.fractions for streamname, streamobject in self.streamregistry.items()]
        return np.array(a)
    
    # Method for outer class "PhysicalProcess"
    #----------------------------------------                
    def _make_numpyarrayof_flowrates(self):
        r""" Creates numpy array of all stream-flowrates for a physical process.
                   
        """
        a = [streamobject.flowrate for streamname, streamobject in self.streamregistry.items()]
        # print("makenumpyFlowrate = ", np.vstack(a))
        return np.vstack(a)

    # Method for outer class "PhysicalProcess"
    #----------------------------------------    
    def _make_numpyarrayof_streamnames(self):
        r""" Creates numpy array of all stream-names for a physical process.
                   
        """
        a = [streamname for streamname, streamobject in self.streamregistry.items()]
        return np.vstack(a)
    
    # Method for outer class "PhysicalProcess"
    #----------------------------------------    
    def _make_arrayof_extrainfo(self):
        r""" Creates numpy array of all stream-extrainfos for a physical process.
                   
        """
        a=[]
        for streamname, streamobject in self.streamregistry.items():
            if hasattr(streamobject, 'extrainfo'):
                a.append(streamobject.extrainfo)
            else:
                a.append([])
        return a
    
    # Method for outer class "PhysicalProcess"
    #-----------------------------------------    
    def _parse_extrainfo(self):
        r"""This method parses the extra info for 'TOTAL flowrate' relationship OR
        'COMPONENT flowrate' relationship between streams.

        Returns
        -------
        parsed_extrainfo : as list of lists
                            [[], [],[]]
            Each sub list contains the following list elements 
          
            CASE 1: COMPONENT flowrate relationship
            [componentindex, left_streamsign, left_streamname, ratio, right_streamname]

            example:
                let extrainfo : ["2:S1=1.2*S4"] 
                This means for component '2', the component flowrate in S1 is 1.2 times that in S4
                So component '2' flowrates in S1 and S4 will be
                    in S1 = F * x (F = total flowrate of S1, and x = component '2' fraction in S1)
                    in S4 = F * x (F = total flowrate of S4, and x = component '2' fraction in S4)
                And equation becomes : abs(F*x for S1) - abs(F*x for S4) = 0
                This equation is NOT formulated by this method/function.
                This method/function ONLY parses the information input by user
                in a manner it can be used later to formulate this equation
                
                So the extrainfo can be parsed as follows:
                componentindex = '2'   
                left_streamsign = +1 if 'S1' is input
                                 -1 if 'S1' is output
                left_streamname = 'S1'
                ratio = '1.2'
                right_streamname = 'S4'
                
                NOTE: "left_streamsign" is NOT needed to formulate the desired
                equation for 'COMPONENT flowrate' relationship. But it is required for
                'OVERALL flowrate' relationship between streams. Therefore to maintain
                uniformity in signature of 'return value', we keep this field

            CASE 2: OVERALL flowrate relationship
            [componentindex, left_streamsign, left_streamname, ratio, right_streamname]
            example:
                let extrainfo : ["S1=1.2*S4"] 
                then:
                componentindex = FALSE
                left_streamsign = +1 if 'S1' is input
                                 -1 if 'S1' is output
                left_streamname = 'S1'
                ratio = 1.2
                right_streamname = 'S4'
                
                Out put of this method is then used to create the equation
                
                left_streamflow - left_streamsign * abs(eval(ratio) * right_streamflow)
                
                The reasoning is that the righthandside should evaluate
                to flowrate of S1, and the sign should be +1 or -1 based on convention
                chosen of input or output streams respectively

        """
        parsed_extrainfo=[]
        for streamname, streamobject in self.streamregistry.items():
            if hasattr(streamobject, 'extrainfo'):
                # print(streamobject.extrainfo)
                for items in streamobject.extrainfo:
                    # first parse for 'component flowrate'
                    # if info is not for 'component flowrate'
                    # the first statement of 'try' will fail
                    # it then helps to decide that the relationship
                    # is for total stream flowrates
                    try:
                        componentindex, other_1 = items.split(":")
                        left_streamname, other_2 = other_1.split("=")          
                        ratio, right_streamname = other_2.split("*")
                        # Check if left_stream and right_stream are both
                        # attached to the same process.
                        # If there are more than one processes that are
                        # interconnected, say process1 and process2,
                        # then the extra information maybe between a stream
                        # in process1 and another stream from process2.
                        # In this case unless a balance is made by placing
                        # a boundary that encloses both these processes
                        # (like an overall balance), then this extra information
                        # is not usable.
                        stream_names_in_process = list(self.streamregistry.keys())
                        if left_streamname in stream_names_in_process and right_streamname in stream_names_in_process:
                            left_streamsign = self._getsignof(self.streamregistry[left_streamname].flowrate)
                            right_streamsign = self._getsignof(self.streamregistry[right_streamname].flowrate)
                            temp = [componentindex, left_streamsign, left_streamname, ratio, right_streamsign, right_streamname]
                            parsed_extrainfo.append(temp)
                    except ValueError:
                        # print(items)
                        left_streamname, other_2 = items.split("=")
                        ratio, right_streamname = other_2.split("*")
                        # Check if left_stream and right_stream are both
                        # attached to the same process.
                        # If there are more than one processes that are
                        # interconnected, say process1 and process2,
                        # then the extra information maybe between a stream
                        # in process1 and another stream from process2.
                        # In this case unless a balance is made by placing
                        # a boundary that encloses both these processes
                        # (like an overall balance), then this extra information
                        # is not usable.
                        stream_names_in_process = list(self.streamregistry.keys())
                        if left_streamname in stream_names_in_process and right_streamname in stream_names_in_process:
                            left_streamsign = self._getsignof(self.streamregistry[left_streamname].flowrate)
                            right_streamsign = self._getsignof(self.streamregistry[right_streamname].flowrate)
                            temp = [False, left_streamsign, left_streamname, ratio, right_streamsign, right_streamname]
                            # print("extra Info", temp)
                            parsed_extrainfo.append(temp)
        return parsed_extrainfo



    # Method for outer class "PhysicalProcess"
    #-----------------------------------------
    
    def _make_numpyarraysof_streamdata(self):
        r""" Creates numpy array of all stream attributes for a physical process.
                   
        """
        # organize information given by user for the process streams
        # into individual arrays for convenient processing of data
        self.npstreamnames = self._make_numpyarrayof_streamnames()
        self.npfractions = self._make_numpyarrayof_fractions()
        self.npflowrates = self._make_numpyarrayof_flowrates()
        self.extrainfos = self._make_arrayof_extrainfo()
        return None

        
    # Method for outer class "PhysicalProcess"
    #----------------------------------------  
    @staticmethod
    def _find_whereplusFs(arg):
        r"""Finds location of "F" in the streams of a physical process.
                   
        """
        # arg : Input parameter = list of flowrates
        a = np.where(arg=="F")
        return a
    
    # Method for outer class "PhysicalProcess"
    #----------------------------------------      
    @staticmethod
    def _find_whereminusFs(arg):
        r"""Finds location of "-F" in the streams of a physical process.
                   
        """
        # arg : Input parameter = list of flowrates
        a = np.where(arg=="-F")
        return a
    
    # Method for outer class "PhysicalProcess"
    #----------------------------------------      
    @staticmethod
    def _find_wherexs(arg):
        r"""Finds location of "x" in the streams of a physical process.
                   
        """
        # arg : Input parameter = list of fractions
        #                         of single stream
        a = np.where(arg=="x")
        return a
   
    # Method for outer class "PhysicalProcess"
    #-----------------------------------------   
    def __str__(self, print_thesestreams=None):
        r"""Implements 'print' function for a physical process.
                  
        """
        if print_thesestreams is None:
            streamstoprint = self.streamregistry
        if print_thesestreams is not None:
            streamstoprint = print_thesestreams            
        # print header info before printing data
        headerstring = "Process streams for " + self.processname.upper() + " are :"
        headerstring_charactercount = (len([chars for items in headerstring for chars in items]))
        print("=" * headerstring_charactercount)
        print(headerstring)
        print("=" * headerstring_charactercount)
        
        # func for mapping fractions of streams for formatting
        # This is done because fractions can be float or string (i.e. 'x')
        def _get_formatspecfor_fractions(arg):
            if type(arg) is str:
                return "{:^6s}".format(arg)
            else:
                return "{:^.4f}".format(arg)
        
        # find max length of string names to set length of display field
        max_streamname_length = max([len(names) for names in streamstoprint.keys()])
        formatfor_streamname = "{0:"+str(max_streamname_length)+"s}"

        # find max length of flow rates to set length of display field
        maxflow1 = max([len(str(round(vals.flowrate,2))) for vals in streamstoprint.values() if type(vals.flowrate) != str], default=0)
        maxflow2 = max([len(vals.flowrate) for vals in streamstoprint.values() if type(vals.flowrate) == str], default=0)
        maxflow = max(maxflow1, maxflow2)
        for streamname, streamobject in streamstoprint.items():
            # prepare flowrate for printing
            # check if flowrate value for a stream is 'F' or '-F' or float
            # and accordingly prepare format specification
            if type(streamobject.flowrate) == str:
                flowformatter = "{1:>"+str(maxflow+3)+"s}" # 3 added to maxflow to match ... see below
            else:
                flowformatter = "{1:>" + str(maxflow+3) + ".2f}" # 3 added to maxflow = one for period, and two for 2 digits after period
            # prepare fractions for printing
            # check if fraction value for a stream is 'x' or float
            # and accordingly prepare format specification
            # use map function to process all fractions
            streamfractions = list(map(_get_formatspecfor_fractions,streamobject.fractions))          
            # prepare full string format
            prtString = "Stream = " + formatfor_streamname + "; Flowrate = " + flowformatter +"; Fractions = {2}; Extra Info = {3} "
            # prepare string values to be printed based on whether 'extrainfo' exists for a stream or not
            if hasattr(streamobject, "extrainfo") :
                print(prtString.format(streamname, streamobject.flowrate, streamfractions, streamobject.extrainfo))
            else:
                print(prtString.format(streamname, streamobject.flowrate, streamfractions, "None"))
        return ("-" * headerstring_charactercount +"\n")*1 + "END\n"
                        # this is symbolic return because __str__ must return something
                        # otherwise print gives TypeError 

        
    # Method for outer class "PhysicalProcess"
    #----------------------------------------    
    def find_unknownfractions(self, streamname=None): # Σx = 1 balance
        r"""Compute unknown fractions 'x'.
        
        Parameters
        ----------
        physical process : Instance of ~.PhysicalProcess
            Physical process for which 'x' are to be computed
        streamname : instance of ~.PhysicalProcess.Stream, default = None
            A specific stream of the PhysicalProcess for which 'x' is to
            be computed. If this parameter is not provided, all streams
            of the physical process are considered for finding the 'x'.
            
                              
        
        Returns
        -------
        A list of lists. Each sublist item can take one of the following form.
        
        - If 'x' is computable : [True, streamname, location_wherexunknown,  missingxvalue]
        - If 'x' are not computable : [False, streamname, textmessage]
        
        
        Notes
        -----
        Attempt is made to compute the unknown 'x' fraction for each stream
        using  Σx = 1. This approach can only work if there is just one
        unknown 'x' in a given stream.
                             
        
        Examples
        --------
        .. code-block:: python

            >>> from pychemengg.massbalances import physicalmassbalance as pmb
            
            # Next create instance of physical process, say it is called 'process1'
            >>> process1 = pmb.PhysicalProcess("process1")

            # Attach three inlet streams to process1.
            >>> process1.S1 = process1.attachstreams(streamnames=["S1"], flowrates=[300], fractions=[[0.1, 0.3, 0.6]])
            >>> process1.S2 = process1.attachstreams(streamnames=["S2"], flowrates=[200], fractions=[[0.1, 0.12, "x"]])
            >>> process1.S3 = process1.attachstreams(streamnames=["S3"], flowrates=[700], fractions=[[0.1, "x", "x"]])
            
            # Apply the find_unknownfractions() method on the process
            >>> unknownx = process1.find_unknownfractions()
            >>> print(unknownx)
            Output is:
            [[False, 'S1', 'There are no missing fractions.'], [True, 'S2', 2, 0.78], [False, 'S3', "More than one 'x' unknown fractions in the stream. Σx = 1 is insufficient to find them."]]
                
            # Apply the find_unknownfractions() method on one of the process stream
            >>> unkn = process1.find_unknownfractions(streamname="S1")
            >>> print(unkn)
            Output is:
            [[False, 'S1', 'There are no missing fractions.']]
  
        
        """
        # Case: Apply to all streams
        unknownx_detailsfound = []
        if streamname == None:
            for strname, streamobject in self.streamregistry.items():
                missingx = streamobject._calc_missing_xfraction()
                if missingx[0] == True:
                    location_wherexunknown = [i for i, val in enumerate(streamobject.fractions) if val =="x"]
                    unknownx_detailsfound.append([True, strname, location_wherexunknown[0], missingx[1]])
                if missingx[0] == False:
                    unknownx_detailsfound.append([False, strname, missingx[1]])

        # Case: Apply balance to single stream when 'streamname' provided in the method argument
        elif streamname != None:
            streamobject = self.streamregistry[streamname]
            missingx = streamobject._calc_missing_xfraction()
            if missingx[0] == True:
                location_wherexunknown = [i for i, val in enumerate(streamobject.fractions) if val =="x"]
                unknownx_detailsfound.append( [True, streamname, location_wherexunknown[0], missingx[1]] )
            if missingx[0] == False:
                unknownx_detailsfound.append( [False, streamname, missingx[1]] )
        return unknownx_detailsfound
          

    # Method for outer class "PhysicalProcess"
    #-----------------------------------------    
    def update_streamfractions(self, new_fractioninfo):
        r"""Update mole/mass fraction of a stream.
        
        This method should be used after finding unknown 'x' using
        'find_unknownfractions'.
        
        
        Parameters
        ----------
        new_fractioninfo : `list`
            A list containing information to update stream fractions.
            
            The form of the list is:
            - [True, streamname, location_wherexunknown, missingx_value]
            - [False, streamname, missingx_value]
        
        
        Returns
        -------
        None
        
       
        Examples
        --------
        .. code-block:: python

            >>> from pychemengg.massbalances import physicalmassbalance as pmb
            
            # Next create instance of physical process, say it is called 'process1'
            >>> process1 = pmb.PhysicalProcess("process1")

            # Attach three inlet streams to process1.
            >>> process1.S1 = process1.attachstreams(streamnames=["S1"], flowrates=[300], fractions=[[0.1, 0.3, 0.6]])
            >>> process1.S2 = process1.attachstreams(streamnames=["S2"], flowrates=[200], fractions=[[0.1, 0.12, "x"]])
            >>> process1.S3 = process1.attachstreams(streamnames=["S3"], flowrates=[700], fractions=[[0.1, "x", "x"]])
            
            # Apply the find_unknownfractions() method on the process
            >>> unknownx = process1.find_unknownfractions()
            >>> print(unknownx)
            Output is:
            [[False, 'S1', 'There are no missing fractions.'], [True, 'S2', 2, 0.78], [False, 'S3', "More than one 'x' unknown fractions in the stream. Σx = 1 is insufficient to find them."]]
            
            # Use unknownx to update streams
            >>> process1.update_streamfractions(unknownx)
            Output is :
            For  " process1 " : the stream " S1 " There are no missing fractions.
            For  " process1 " : the stream " S2 " The new component fraction @ position = " 2 " is now =  0.78
            For  " process1 " : the stream " S3 " More than one 'x' unknown fractions in the stream. Σx = 1 is insufficient to find them.
            
        
        
        See Also
        --------
        pychemengg.massbalances.physicalmassbalance.PhysicalProcess.find_unknownfractions
        
        """
        for info in new_fractioninfo:
            if info[0] == True:
                self.streamregistry[info[1]].fractions[info[2]]=info[3]
                print("For ", '"', self.processname, '"', ": the stream", '"', info[1], '"', "The new component fraction @ position =", '"',info[2],'"', "is now = ", info[3])
                
            if info[0] == False:
                print("For ", '"', self.processname, '"', ": the stream", '"', info[1], '"', info[2])
 
    # Method for outer class "PhysicalProcess"
    #----------------------------------------- 
    def find_unknownflowrate(self): # ΣFin + ΣFout = 0
        r"""Compute unknown flowrate "F" or "-F".
        
        Parameters
        ----------
        physical process : Instance of ~.PhysicalProcess
            Physical process for which 'F'/'-F' are to be computed.            
                              
        
        Returns
        -------
        A list. Each list can take one of the following form.
        
        - If 'F'/'-F' is computable : [[True, streamname, missingflowrate_value]]
        - If 'F'/'-F' is not computable : [[False, textmessage]]
        
        
        Notes
        -----
        Attempt is made to compute the unknown flowrate for a process
        using  # ΣFin + ΣFout = 0. This approach can only work if there is
        just one unknown 'F' or '-F' in a given process.
                             
        
        Examples
        --------
        .. code-block:: python

            >>> from pychemengg.massbalances import physicalmassbalance as pmb
            
            # Next create instance of physical process, say it is called 'process1'
            >>> process1 = pmb.PhysicalProcess("process1")

            # Attach two inlet streams and one outlet stream to process1.
            >>> process1.S1 = process1.attachstreams(streamnames=["S1"], flowrates=[300], fractions=[[0.1, 0.3, 0.6]])
            >>> process1.S2 = process1.attachstreams(streamnames=["S2"], flowrates=[200], fractions=[[0.1, 0.12, "x"]])
            # The stream below is the outlet stream. Note the flowrate has a negative sign.
            >>> process1.S3 = process1.attachstreams(streamnames=["S3"], flowrates=["-F"], fractions=[[0.1, "x", "x"]])
            
            # Apply the find_unknownflowrate() method on the process
            >>> unknownF = process1.find_unknownflowrate()
            >>> print(unknownF)
            Output is:
            [[True, 'S3', -500.0]]
                        
        """
        unknownflowrate_detailsfound = []
        plusF_location, minusF_location, x_location = self._find_unknownxandFlocations()
        unknownplusF_count = len(plusF_location[0])
        unknownminusF_count = len(minusF_location[0])
        if (unknownplusF_count + unknownminusF_count) == 1:
            for flow in self.npflowrates:
                knownflowrates = [float(flow[0]) for flow in self.npflowrates if (flow != "F" and flow != "-F")]
                sumof_knownflowrates = np.sum(np.array(knownflowrates).astype(float))
                # Fsum = np.sum(np.array([float(flow[0]) for flow in self.npflowrates if (flow != "F" and flow != "-F")]).astype(float))
                missingflowrate = (0 - sumof_knownflowrates) # because ΣF = 0
                if unknownplusF_count == 1:
                    streamwith_unknownflowrate = self.npstreamnames[(plusF_location)]
                if unknownminusF_count == 1:
                    streamwith_unknownflowrate = self.npstreamnames[(minusF_location)]
                # print("name and missing flowrate", streamwith_unknownflowrate, missingflowrate) 
                unknownflowrate_detailsfound.append([True, streamwith_unknownflowrate[0], missingflowrate])
                return unknownflowrate_detailsfound
        
        if (unknownplusF_count + unknownminusF_count) > 1:
            text = "More than one 'F'/'-F' unknown flowrates in stream. Try the 'solvesystem' method"
            unknownflowrate_detailsfound.append([False, text])
            return unknownflowrate_detailsfound
        
        if (unknownplusF_count + unknownminusF_count) == 0.0:
            text = "There are no missing flowrates"
            unknownflowrate_detailsfound.append([False, text])
            return unknownflowrate_detailsfound


    # Method for outer class "PhysicalProcess"
    #--------------------------------------------     
    def update_streamflowrates(self, newflow):
        r"""Update flowrate of a stream.
        
        This method should be used after finding unknown 'x' using
        'find_unknownflowrate'.
        
        
        Parameters
        ----------
        new_flow : `list`
            A list containing information to update stream fractions.
            
            The form of the list is:
            - [[False, textmessage]]
            - [[True, streamname, missingflowrate_value]]
                              
        
        Returns
        -------
        None
        
       
        Examples
        --------
        .. code-block:: python

            >>> from pychemengg.massbalances import physicalmassbalance as pmb
            
            # Next create instance of physical process, say it is called 'process1'
            >>> process1 = pmb.PhysicalProcess("process1")

            # Attach two inlet streams and one outlet stream to process1.
            >>> process1.S1 = process1.attachstreams(streamnames=["S1"], flowrates=[300], fractions=[[0.1, 0.3, 0.6]])
            >>> process1.S2 = process1.attachstreams(streamnames=["S2"], flowrates=[200], fractions=[[0.1, 0.12, "x"]])
            # The stream below is the outlet stream. Note the flowrate has a negative sign.
            >>> process1.S3 = process1.attachstreams(streamnames=["S3"], flowrates=["-F"], fractions=[[0.1, "x", "x"]])
            
            # Apply the find_unknownflowrate() method on the process
            >>> unknownF = process1.find_unknownflowrate()
            >>> process1.update_streamflowrates(unknownF)
            Output is:
            For  " process1 " : the stream " S3 " has new flowrate = -500.0
        
        
        See Also
        --------
        pychemengg.massbalances.physicalmassbalance.PhysicalProcess.find_unknownflowrate
        
        """
        # print(newflow)
        for streamdata in newflow:
            if streamdata[0] == True:
                self.streamregistry[streamdata[1]].flowrate=streamdata[2]
                print("For ", '"', self.processname, '"', ": the stream",'"', streamdata[1], '"', "has new flowrate =", streamdata[2])

            if streamdata[0] == False:
                print("For ", '"', self.processname, '"', streamdata[1])
                


    # Method for outer class "PhysicalProcess"
    #-----------------------------------------
    def perform_component_massbalance(self):
        r"""Execute component mass balance on a process.
        
        :math:`(ΣF_ix_{i,component})` = 0; where i = stream index
        
       
        Parameters
        ----------
        physical process : Instance of ~.PhysicalProcess
            Physical process for which component mass balance is to be performed.       
                              
        
        Returns
        -------
        None
        
        
        See Also
        --------
        pychemengg.massbalances.physicalmassbalance.PhysicalProcess.find_unknownflowrate
        pychemengg.massbalances.physicalmassbalance.PhysicalProcess.find_unknownfractions
        
        
        
        Notes
        -----
        **Logic for component mass balance**

        
        The nparray of stream flowrates has shape (n,1):
            
        - n : number of streams
                

        The nparray of stream flowrates has shape (n,m): 
            
        - n : number of streams
        - m : number of components per stream
                

        
        For example consider three streams:
        
        ==================   ============   =================
        Stream               Flowrate       Fractions
        ==================   ============   =================
        Inlet 1              F1             [x11, x12, x13]
        Inlet 2              F2             [x21, x22, x23]
        Outlet 3             F3             [x31, x32, x33]
        ==================   ============   =================

        The nparray for flowrates is:
            
            A = 
            
            [[F1]
             
            [F2]
             
            [F3]] 
             
            shape = (3, 1)
        
        
        The nparray for component fractions is:
            
            B =
            
            [[x11, x12, x13]
            
            [x21, x22, x23]
            
            [x31, x32, x33]] 
            
            shape = (3, 3)
        
        
        :math:`x_{m,n}` : 
            
            - m : stream index
            - n : component index
        
        To perform component mass balance, the row of matrix 'A' is multiplied
        to row of matrix 'B' 
        
        C = 
        
        [[F1*x11, F1*x12, F1*x13]
        
        [F2*x21, F2*x22, F2*x23]
        
        [F3*x31, F3*x32, F3*x33]]
        
        Each column of 'C' represents flowrate of a component.
        By summing the column of 'C' the component mass balance can be
        performed. 
        
        - Sum of column 1 = Flowrate balance for component 1
        - Sum of column 2 = Flowrate balance for component 2
        
        Examples
        --------
        .. code-block:: python

            >>> from pychemengg.massbalances import physicalmassbalance as pmb
            
            # Next create instance of physical process, say it is called 'process1'
            >>> process1 = pmb.PhysicalProcess("process1")

            # Attach two inlet streams and one outlet stream to process1.
            >>> process1.S1 = process1.attachstreams(streamnames=["S1"], flowrates=[300], fractions=[[0.4333, 0.3, 0.2667]])
            >>> process1.S2 = process1.attachstreams(streamnames=["S2"], flowrates=[200], fractions=[[0.1, 0.12, 0.78]])
            # The stream below is the outlet stream. Note the flowrate has a negative sign.
            >>> process1.S3 = process1.attachstreams(streamnames=["S3"], flowrates=[-500], fractions=[[0.3, 0.228, 0.472]])
            
            # Apply the perform_component_massbalance() method on the process
            >>> componentbalance = process1.perform_component_massbalance()
            >>> print("component balance: ", componentbalance)
            Output is:
            component balance:  [-0.01  0.    0.01]
            # This means that:
            # For component 1: ΣFx = -0.01
            # For component 2: ΣFx = 0.
            # For component 3: ΣFx = 0.01
            # The above balance can be set up as part of a function to solve
            # for unknown "F", "-F", and "x". See 'User Tutorial'
        
        """
        # Performs (ΣFx )in = (ΣFx)out for each component
        # Returns numpy array of mass balance for each component
        # if the process is uniquely defined then these terms will
        # be close to zero
        self._make_numpyarraysof_streamdata()
        # print(self.npflowrates, "\n shape =", self.npflowrates.shape)
        # print("\n",self.npfractions, "\n shape =", self.npfractions.shape)
        npcomponent_massbalance = (self.npflowrates * self.npfractions).sum(axis=0)
        
        # THIS IS ALTERNATE CODE
        #------------------------
        # npstream_flowandfraction_data = [[np.array(streamobject.flowrate),np.array(streamobject.fractions)]
        #                for streamname, streamobject in self.streamregistry.items()]
        # component_mass = [numpyFlow * numpyFractions for numpyFlow, numpyFractions in npstream_flowandfraction_data]
        # npcomponent_mass = np.array(component_mass)
        # npcomponent_massbalance = npcomponent_mass.sum(axis=0) 
        
        return npcomponent_massbalance # return type <class 'numpy.ndarray'> ... 1 x n array
    
    
    # Method for outer class "PhysicalProcess"
    #----------------------------------------- 
    def _find_unknownxandFlocations(self):
        r"""Find where 'F', '-F', 'x' are located in different streams.
        
        Parameters
        ----------
        None
        

        Returns
        -------
        plusF_location : tuple (nparray of `int`)
            Column numbers of array 'npflowrates' where "F" is located.
        minusF_location : tuple (nparray of `int`)
            Column numbers of array 'npflowrates' where "-F" is located.
        x_location : tuple (nparray of `int`, nparray of `int`)
            First nparray of tuple has row index, second nparray of tuple has column index

        """
        self._make_numpyarraysof_streamdata()
        plusF_location = self._find_whereplusFs(self.npflowrates)
        minusF_location = self._find_whereminusFs(self.npflowrates)
        x_location = self._find_wherexs(self.npfractions)
        return (plusF_location, minusF_location, x_location)  
    
    
    # Method for outer class "PhysicalProcess"
    #-----------------------------------------   
    def _getsignof(self, F):
        r"""Find sign of the flowrate of a stream.
        
        """
        if isinstance(F, str):
            if F == "F":
                return 1
            if F == "-F":
                return -1
        if (isinstance(F, int) or isinstance(F, float)):
            if F > 1:
                return 1
            if F < 1:
                return -1
            
            
    # Method for outer class "PhysicalProcess"
    #-----------------------------------------  
    def degreesoffreedom(self):
        r"""Perform degree of freedom analysis on the physical process.
        
        Parameters
        ----------
        None : Method is applied to the physical process.
        
        Returns
        -------
        None : None. 
            Multiple "self" attributes are created. These are used to
            formulate balances that are required to solve for "F", "-F", and "x"
        
        Attributes
        ----------
        
        unknown_x_count : `int`
            Count of total unknown fractions in streams attached to the process.
        unknownplusF_count : `int`
            Count of total unknown inlet stream flowrates.
        unknownminusF_count : `int`
            Count of total unknown outlet stream flowrates.
        unknown_Fs_count : `int`
            Count of total unknown inlet + outlet stream flowrates.
        totalunknowns : `int`
            Sum of unknown_Fs_count + unknown_x_count
        extrainfo_count : `int`
            Count of total extra information available for the process.
        component_count : `int`
            Total components in the process streams
        streamswith_unknownx_count : `int`
            Count of streams that contain at least one unknown component fraction.
        totalequations_possible : `int`
            Sum of extrainfo_count + component_count + streamswith_unknownx_count
      
        """
        self._make_numpyarraysof_streamdata()
        self.plusF_location, self.minusF_location, self.x_location = self._find_unknownxandFlocations()
        self.unknownx_count = len(self.x_location[0])
        self.unknownplusF_count = len(self.plusF_location[0])
        self.unknownminusF_count = len(self.minusF_location[0])
        self.unknown_x_count = self.unknownx_count
        self.unknown_Fs_count = self.unknownplusF_count + self.unknownminusF_count
        self.totalunknowns = self.unknown_x_count + self.unknown_Fs_count
        self.extrainfo_records = self._parse_extrainfo() 
        self.extrainfo_count  = len(self.extrainfo_records)
        self.streamswith_unknownx_count = len(set(self.x_location[0]))
        self.component_count = len(self.npfractions[0])        
        # When there are multiple processes, it is possible
        # some of them may have fewer components in their streams.
        # Find if a component(s) has "0" fractions in all streams.
        # If, yes, subtract it from self.component_count
        iscolumn_with_all_xfractions_equalto_0 = np.all((self.npfractions == '0'), axis=0)
        countofcolumns_with_all_xfractions_equalto_0  = np.count_nonzero(iscolumn_with_all_xfractions_equalto_0 == True)
        if countofcolumns_with_all_xfractions_equalto_0 > 0:
            self.component_count = self.component_count - countofcolumns_with_all_xfractions_equalto_0
        self.totalequations_possible = self.extrainfo_count + self.component_count + self.streamswith_unknownx_count

        beginning_text = "Degrees of freedom analysis for " + self.processname.upper()
        print("=" * (len(beginning_text)))
        print(beginning_text)
        print("=" * (len(beginning_text)))
        print("Number of unknown flowrates :-->", self.unknown_Fs_count)
        print("Number of unknown 'x' fractions :-->", self.unknown_x_count)
        print("Total unknowns :-->", self.totalunknowns)
        print("-"*20)
        print("Total possible component balances (ΣFx)in = (ΣFx)out :--> ", self.component_count)
        print("Total possible sum of stream fractions is unity balances (Σx) = 1 :-->", self.streamswith_unknownx_count)
        print("Other extra equations :-->", self.extrainfo_count)
        print("Total equations :-->", self.totalequations_possible)

        if self.totalunknowns > self.totalequations_possible:
            print("\nSystem is underspecified because there are more unknowns than available equations")
            print("\nThere are = " , self.totalunknowns , "unknowns and" , self.totalequations_possible , "possible equations")
            ending_text = "End of degree of freedom analysis for " + self.processname.upper()
            print(ending_text)
            print("-" * len(ending_text))
            print() # leave an empty line
            return None
        
        if self.totalunknowns == self.totalequations_possible:
            print("\nSystem can be solved")
            print("\nThere are = " , self.totalunknowns , "unknowns and" , self.totalequations_possible , "possible equations")
            ending_text = "End of degree of freedom analysis for " + self.processname.upper()
            print("-" * len(ending_text))
            print(ending_text)
            print("-" * len(ending_text))
            print("\n") # leave an empty line
            return None
        
        if self.totalunknowns < self.totalequations_possible:
            print("\nIt maybe possible to solve the system")
            print("\nThere are = " , self.totalunknowns , "unknowns and" , self.totalequations_possible , "possible equations")
            ending_text = "End of degree of freedom analysis for " + self.processname.upper()
            print("-" * len(ending_text))
            print(ending_text)
            print("-" * len(ending_text))
            print("\n") # leave an empty line
            return None
    
    
    # Method for outer class "PhysicalProcess"
    #-----------------------------------------
    @staticmethod
    def _create_solvedsystem_summary(processname, npstreamnames,
                                     npflowrates, npfractions, extrainfos):
        r"""Creates the solution object of the physical process being solved.
        
        The object created is an instance of PhysicalProcess class.
        Different streams tat are originally associated with the process
        that is being solved get attached to this object with all known 
        quantities.
        """
        string = "SOLUTION TO "  + processname.upper()
        solvedsystem = PhysicalProcess(string)
        for num, strmName in enumerate(npstreamnames):
            
            # The following is a hack to attachstreams and create .streamname attribute
            #--------------------------------------------------------------------------
            # headerstring = "solvedsystem."+strmName[0]+"= solvedsystem.attachstreams(streamnames=["+"""strmName[0]]"""+" )"
            # exec(headerstring)
            # Then Assign flowrate and fractions to streamobject
            # solvedsystem.streamregistry[strmName[0]].flowrate = npflowrates[num][0]
            # solvedsystem.streamregistry[strmName[0]].fractions = npfractions[num]
            
            # if len(extrainfos[num]) > 0:
            #     solvedsystem.streamregistry[strmName[0]].extrainfo = extrainfos[num]
            
            # Following is a much cleaner way
            #----------------------------------
            # First create streams
            streamobject = solvedsystem.attachstreams(streamnames=[strmName[0]], 
                                                      flowrates=npflowrates[num], 
                                                      fractions=[npfractions[num]],
                                                      extrainfos=[extrainfos[num]])
            
            # Reason why [  ] is placed around some parameters
            #--------------------------------------------------            
            # print(strmName[0]) :--> output is 'S2' 
            # streamnames needs a LIST as input, so use [ strmName[0] ]
            
            # print(npflowrates[num]) :--> output is [-3972.33]
            # flowrates needs a LIST, and this is already a list so NO CHANGE
            
            # print(npfractions[num]) :--> output is [0.0, 1.0]
            # fractions needs a LIST of LIST as input, so use [ npfractions[num] ]
            
            # print(extrainfos[num]) :--> output is []
            # extrainfos needs a LIST of LIST as input, so use [ extrainfos[num] ]
            
            # Next set attribute such as mixer1.S1 = streamobject
            setattr(solvedsystem, strmName[0], streamobject)
        return solvedsystem
        

    # Method for outer class "PhysicalProcess"
    #-----------------------------------------           
    def solvesystem(self):
        r"""Solves the physical process to find "F", "-F", and "x"
        
        
        Parameters
        ----------
        None : Apply the method to the physical process to be solved.
        

        Returns
        -------
        solution : Instance of the class PhysicalProcess
            A new instance of PhysicalProcess is created that
            stores the solution. The solution can be printed.
            
        """
        self.degreesoffreedom() 
        # Perform degrees of freedom analysis
        # This executes code that will create
        # many of the "self.variables" that are needed for
        # analysis below.

        # Method inner to solvesystem()
        #-------------------------------
        def _assemble_all_possible_balance_and_extrainfo_equations(guess):
            # This function is called to assemble and return
            # all possible balance equations and extrainfo equations.
            # If the degrees of freedom is zero, then all of them
            # are used to solve the system. If the degrees of freedom
            # is such that number of equations are greater than number
            # of unknowns then the system is overspecified, however,
            # still, if the equations are consistent, they can be used
            # to compute the unknowns. For this, it is important to
            # weed out equations that express the same relationship
            # and are therfore NOT UNIQUE. Only one of these "non-unique"
            # equations can be used. Furthermore, TRIVIAL equations
            # for which F*x terms in a component mass balance vanish
            # if x=0 for all the unknown "F" or "-F" terms,
            # should be weeded out.
            # This function only assembles all possible equations
            # and returns them to another function where this weeding out
            # is performed.
            
            # If system is underspecified, print message ---> EXIT.
            if self.totalunknowns > self.totalequations_possible:
                print("\nSystem is underspecified because there are more unknowns than available equations")
                print("\nThere are = ", self.totalunknowns, "unknowns and", self.totalequations_possible,
                      "equations")
                return None
            
            # If total unknowns <= total equations ---> PROCEED AHEAD
            if self.totalunknowns <= self.totalequations_possible:
                guessindex = 0 # initialize index to array of guess value
                
                # (A) ASSIGN GUESS VALUES TO UNKNOWN VARIABLES
                # ---------------------------------------------
                # 1) First assign 'guess' values to unknown 'x' 
                if self.unknownx_count > 0:
                    self.npfractions[(self.x_location)] = guess[guessindex : self.unknownx_count]
                    guessindex = guessindex + self.unknownx_count
                self.npfractions = self.npfractions.astype(float)
                # 2) Next assign 'guess' values to unknown 'F'
                if self.unknownplusF_count > 0:
                    self.npflowrates[(self.plusF_location)] = abs(guess[guessindex : guessindex + self.unknownplusF_count])
                    # NOTE: absolute value is used on RHS because "F"
                    # is supposed to be positive, and if the guess is
                    # negative, the solution may not converge.
                    guessindex = guessindex + self.unknownplusF_count
                # 3) Next assign 'guess' values to unknown '-F'
                if self.unknownminusF_count > 0:
                    self.npflowrates[(self.minusF_location)] = -abs(guess[guessindex : guessindex + self.unknownminusF_count])
                    # NOTE:  "-" absolute value is used on RHS because "F"
                    # is supposed to be negative, and if the guess is
                    # positive, the solution may not converge.
                    guessindex = guessindex + self.unknownminusF_count
                self.npflowrates = self.npflowrates.astype(float)
                
                # (B) NEXT CREATE EQUATIONS TO ALLOW SOLUTION OF 'x' and 'F'/'-F'
                # --------------------------------------------------------------------
                
                # ----------------------------------------------                
                # # 1) Start by forming equations using'Extra Info' of streams
                # ----------------------------------------------
                # if self.extrainfo_count > 0:
                #     extrainfoindex = 0 # initialize index into extrainfo relations
                #     extrarelations ={} # to store extra relation equations
                #     # Process extra info records
                #     for indx, relation in enumerate(self.extrainfo_records):
                #         isextrainfo_between_component_flowrates = relation[0]
                #         left_streamsign = relation[1]
                #         left_streamname = relation[2]
                #         ratio = relation[3]
                #         right_streamsign = relation[4]
                #         right_streamname = relation[5]
                #         left_streamflow = self.npflowrates[np.where(self.npstreamnames == left_streamname)]
                #         right_streamflow = self.npflowrates[np.where(self.npstreamnames == right_streamname)]
                #         if isextrainfo_between_component_flowrates == False: # then extrainfo is between components of streams
                #             if left_streamsign * right_streamsign > 0: # both streams are either '+':inlet or '-':outlet
                #                 extrarelations[extrainfoindex] = left_streamsign*abs(left_streamflow) - right_streamsign * abs(((eval(ratio) * right_streamflow)))
                #                 extrainfoindex = extrainfoindex + 1
                #             if left_streamsign * right_streamsign < 0: # one stream is '+':inlet other is '-':outlet
                #                 extrarelations[extrainfoindex] = (left_streamsign*(left_streamflow)) - abs(right_streamsign * (((eval(ratio) * right_streamflow))))
                #                 extrainfoindex = extrainfoindex + 1
                #         if isextrainfo_between_component_flowrates != False: # then extra relationship is between total stream flowrates
                #             column_coordinate = int(isextrainfo_between_component_flowrates[0])-1
                #             left_rowcoordinate = (np.where(self.npstreamnames == left_streamname))
                #             right_rowcoordinate = (np.where(self.npstreamnames == right_streamname)) 
                #             left_streamfraction = self.npfractions[left_rowcoordinate[0][0], column_coordinate]
                #             right_streamfraction = self.npfractions[right_rowcoordinate[0][0], column_coordinate]
                #             if left_streamsign * right_streamsign > 0: # both streams are either '+':inlet or '-':outlet
                #                 extrarelations[extrainfoindex] = left_streamsign*abs(left_streamfraction * left_streamflow) - right_streamsign * abs(eval(ratio) * right_streamfraction * right_streamflow) 
                #                 extrainfoindex = extrainfoindex + 1
                #             if left_streamsign * right_streamsign < 0: # one stream is '+':inlet other is '-':outlet
                #                 extrarelations[extrainfoindex] = left_streamsign*(left_streamfraction * left_streamflow) - abs(right_streamsign * eval(ratio) * right_streamfraction * right_streamflow) 
                #                 extrainfoindex = extrainfoindex + 1
                # --------------------------------------------------
                # The code above ensures that equations can be solved.
                # The < 0 and < 0 conditions are imposed
                # However by forcing '-' or '+' sign to guess variables while
                # they are assigned to "-F" or "F" unknowns, this appears
                # to be no longer required. However, until more testing is
                # done, the above code is retained in case it is needed again.
                
                # 1) Start by forming equations using'Extra Info' of streams
                if self.extrainfo_count > 0:
                    extrainfoindex = 0 # initialize index into extrainfo relations
                    extrarelations ={} # to store extra relation equations
                    # Process extra info records
                    for indx, relation in enumerate(self.extrainfo_records):
                        isextrainfo_between_component_flowrates = relation[0]
                        left_streamsign = relation[1]
                        left_streamname = relation[2]
                        ratio = relation[3]
                        right_streamsign = relation[4]
                        right_streamname = relation[5]
                        left_streamflow = self.npflowrates[np.where(self.npstreamnames == left_streamname)]
                        right_streamflow = self.npflowrates[np.where(self.npstreamnames == right_streamname)]
                        if isextrainfo_between_component_flowrates == False: # then extrainfo is between total stream flowrates
                            extrarelations[extrainfoindex] = abs(left_streamflow) - abs(((eval(ratio) * right_streamflow)))
                            extrainfoindex = extrainfoindex + 1                            
                        if isextrainfo_between_component_flowrates != False: # then extra relationship is between flowrates of components
                            column_coordinate = int(isextrainfo_between_component_flowrates[0])-1
                            left_rowcoordinate = (np.where(self.npstreamnames == left_streamname))
                            right_rowcoordinate = (np.where(self.npstreamnames == right_streamname)) 
                            left_streamfraction = self.npfractions[left_rowcoordinate[0][0], column_coordinate]
                            right_streamfraction = self.npfractions[right_rowcoordinate[0][0], column_coordinate]
                            extrarelations[extrainfoindex] = abs(left_streamfraction * left_streamflow) - abs(eval(ratio) * right_streamfraction * right_streamflow) 
                            extrainfoindex = extrainfoindex + 1
                # 2) Create the component/species balance: (ΣFx)in - (ΣFx)out = 0
                componentmassbalance = (self.npflowrates * self.npfractions).sum(axis=0)
                # 3) Create the component fraction balance: (Σx) = 0
                componentmassfractionbalance = 1-(self.npfractions.sum(axis=1))
                
                # (C) NEXT ASSEMBLE EQUATIONS
                # ---------------------------
                # These are equations that contain 'x' and '+F'/'-F' and 
                # will be solved using a solver
                vareqindex = 0 # index for number of equations assembled
                self.var_equations_registry = {} # Registry to store equation information
                # KEY: tuple of :
                # ("extrainfo",count) --> tracks extrainfo equations
                # ("componentbalance", column number) --> tracks component balance equations
                # ("fractionbalance", row number) --> tracks fraction balance equations
                
                # 1) Start by assembling 'EXTRA INFO' equations      
                if self.extrainfo_count > 0:
                    for equation in extrarelations.values():
                        dict_key = ("extrainfo", vareqindex)
                        self.var_equations_registry[dict_key] = equation[0]
                        vareqindex = vareqindex + 1
                # 2) Next assemble equations involving component mass/mole balances
                # ΣFixi = 0 for each component i
                # Assemble this equation ONLY if at least one 'x' or 'F'/'-F'
                # is unknown
                self.columns_used_in_componentbalance_equations = [] # to track columns used
                if self.unknownx_count > 0:
                    rowcontaining_unknownx = self.x_location[0]
                    columncontaining_unknownx = self.x_location[1]
                    # x_location stores 'row' and 'column' numbers where 'x' exists
                    rownumber_withunknownxcount =  Counter(rowcontaining_unknownx)
                    columnnumber_withunknownxcount = Counter(columncontaining_unknownx)
                    x_frequency_column_pair = [(column, frequency) for column, frequency in columnnumber_withunknownxcount.items()]
                    sorted_x_frequency_column_pair = sorted(x_frequency_column_pair, key=lambda frequency: frequency[1], reverse=True)
                    x_frequency_row_pair = [(row, frequency) for row, frequency in rownumber_withunknownxcount.items()]
                    sorted_x_frequency_row_pair = sorted(x_frequency_row_pair, key=lambda frequency: frequency[1], reverse=True)                      
                    for column_location, x_count in sorted_x_frequency_column_pair:
                        if column_location not in self.columns_used_in_componentbalance_equations:
                            dict_key = ("componentbalance", column_location)
                            self.var_equations_registry[dict_key] = componentmassbalance[column_location]
                            self.columns_used_in_componentbalance_equations.append(column_location)
                            vareqindex = vareqindex + 1                
                if self.unknown_Fs_count > 0:    
                    for column, val in enumerate(componentmassbalance):
                        if column not in self.columns_used_in_componentbalance_equations:
                            dict_key = ("componentbalance", column)
                            self.var_equations_registry[dict_key] = val                            
                            self.columns_used_in_componentbalance_equations.append(column)
                            vareqindex = vareqindex + 1
                # 3) Next assemble component FRACTION balances Σx= 1
                # but don't use streams for which all 'x' are known
                self.rows_used_in_fractionbalance_equations=[] # to track rows used
                if self.unknownx_count > 0:
                    for row_location, x_count, in sorted_x_frequency_row_pair:
                        if row_location not in self.rows_used_in_fractionbalance_equations:
                            dict_key = ("fractionbalance", row_location)
                            self.var_equations_registry[dict_key] = componentmassfractionbalance[row_location]
                            self.rows_used_in_fractionbalance_equations.append(row_location)
                            vareqindex = vareqindex + 1
            return self.var_equations_registry
        
        
        # Method inner to solvesystem()
        #-------------------------------
        def _get_column_locations_for_trivial_equations():
            # Check for TRIVIAL equations
            # If 'F' or '-F' is an unknown in a column vector of flowrates
            # then if all the correspoding 'x' in a certain column are zero
            # then componentbalance term F * x = 0, and it vanishes.
            # The resulting equation of component balance is a scalar quantity that is
            # independent of 'F'/'-F' and 'x' and should be excluded from use
            # by solver.
            # Prepare rows where flowrate is "F" or "-F"
            row_withPlusF_or_minusF = set(self.minusF_location[0]).union(set(self.plusF_location[0]))
            # Initialize a placeholder
            iscolumn_xfraction_zero = []   
            # Loop across the row where 'F' or '-F' exists and check
            # if x-fraction is zero or not.
            for row in row_withPlusF_or_minusF:
                g = [True if x == '0' else False for x in self.npfractions[row]]
                iscolumn_xfraction_zero.append(g)
            # If x-fraction at any location is zero, then check if
            # x-fraction is zero in that column for all places where 'F' or '-F'
            # occurs. If yes, then F * x = 0 for any F and equation if trivial
            if len(iscolumn_xfraction_zero) != 0:
                iscolumn_xfraction_zero = np.asarray(iscolumn_xfraction_zero)    
                trivial_equation_column_location = np.all((iscolumn_xfraction_zero==True), axis=0)
                trivial_equation_column_location = [column_number for column_number, boolean in enumerate(trivial_equation_column_location) if boolean == True ]
            else:
                trivial_equation_column_location = []
            return trivial_equation_column_location
            
            

        # Method inner to solvesystem()
        #-------------------------------
        def _get_dict_keys_for_nonunique_equations():
            # This function assigns arbitrary but unique values to 'x', 'F', and '-F'
            # Values for the different equations are computed.
            # These values are compared to each other to determine
            # whether any of them are equal to each other.
            # If they are, this will imply that those two or more equations
            # are not unique and as such only one of them will be used.
            # Create a variable to store values that correspond
            # to 'x', 'F', and '-F'
            values_to_check_equation_uniqueness = np.ones(self.totalunknowns)            
            # Generate unique values to be assigned to for 'x'
            x_checkvalues = np.linspace(0.01, 1, self.unknown_x_count)
            # Generate unique values to be assigned to for 'F' and '-F'
            F_checkvalues = np.linspace(3, 1000, self.unknown_Fs_count)
            # Assign these values
            values_to_check_equation_uniqueness[0:self.unknown_x_count] = x_checkvalues
            values_to_check_equation_uniqueness[self.unknown_x_count:] = F_checkvalues
            # Pass these 'x', 'F' and '-F' values to the equations assembled
            # in function called _assemble_all_possible_balance_and_extrainfo_equations
            all_possible_equations = _assemble_all_possible_balance_and_extrainfo_equations(values_to_check_equation_uniqueness)
            # Analyse the numerical values obtained from each equation.
            # Identify any equation values that are repeated.
            # Repetition of equation result, despite assigning each 'x', 'F' and '-F'
            # a unique value means that those equations must be identical
            # all_possible_equations = {"1":1, "2":1, "3":2, "4":4, "5":3, "11":123, "34":1, "ss":123}
            # create dict of equation values and their frequencies
            equationvalue_frequencies= Counter(all_possible_equations.values())
            # if equations are not unique, value of equation will be repeated
            non_unique_values = [item for item, count in equationvalue_frequencies.items() if count > 1]
            # find dict key corresponding to repeating equation value
            non_unique_keys = [[key for key, value in all_possible_equations.items() if value==itemtomatch] for itemtomatch in non_unique_values]
            # Of these non-unique values, the very first one will be used
            # and the rest will be removed/deleted, so start at index '1' not '0'
            non_unique_keys_todelete = [item[1:] for item in non_unique_keys]
            non_unique_keys_todelete_as_flatlist = [value for item in non_unique_keys_todelete for value in item]
            return non_unique_keys_todelete_as_flatlist
            
        
        def _unique_equations_to_solve(guessvalues, keys_to_delete):
            # This function is called by 'leastsq'
            # The guessvalues are passed to
            # "_assemble_all_possible_balance_and_extrainfo_equations"
            # This creates and returns all possible equations
            all_possible_equations = _assemble_all_possible_balance_and_extrainfo_equations(guessvalues)
            # From these possible equations, the unique ones are selected
            # and returned
            for keys in keys_to_delete:
                del all_possible_equations[keys]
            return list(all_possible_equations.values())[0:self.totalunknowns]
        
        
        def _prepare_guess_values():
            # PREPARE GUESS VALUES TO CALL 'leastsquare' ----> called 'leastsq'
            # -----------------------------------------------------------------
            self.guessvalues = np.ones(self.totalunknowns)
            guessindex = 0
            # guessvalues are used to assign values to unknown 'x' and 'F'
            # First the 'x' values are assigned, and then 'F' values
            # So we construct guessvalues such that it's initial elements
            # corresponding to 'x' are initialized to 0.5 
            # since 'x' lies between 0 and 1
            if self.unknown_x_count > 0:
                for i in range(self.unknown_x_count):
                    self.guessvalues[guessindex] = 0.5
                    guessindex = guessindex + 1
            # The remainder elements of guessvalues must be set equal to flow rates
            # that are already known.
            # So we find known flow rate values to use them in initilization.
            # If we randomly initialize guess values of flowrates to say 1, 100, or 1000 etc
            # then the fsolve function maynot converge because these random values
            # maybe far from the true solution
            knownflowrates = [float(flow[0]) for flow in self.npflowrates if (flow[0] != "F" and flow[0] != "-F")]
            sorted_knownflowrates = sorted(knownflowrates, reverse=True)
            if self.unknown_Fs_count > 0:
                for i in range(self.unknown_Fs_count):
                    if i < len(knownflowrates):
                        self.guessvalues[guessindex] = sorted_knownflowrates[i]
                        guessindex = guessindex + 1
                    else: # use the last flowrate in list as guess value
                        self.guessvalues[guessindex] = knownflowrates[-1]
                        guessindex = guessindex + 1
            return self.guessvalues
        
        # STEPS to Solve System
        # Step 1) Prepare guess values
        _prepare_guess_values()
        # Step 2) Obtain locations of trivial equations
        trivial_equation_column_locations = _get_column_locations_for_trivial_equations()
        # create keys corresponding to the trivial equation columns
        # the key will be of signature ---> tuple ("componentbalance", column)
        keys_for_trivial_equation_column_locations = [("componentbalance", column) for column in trivial_equation_column_locations]        
        # Step 3) Obtain locations for unique equations
        keys_of_nonunique_equations = _get_dict_keys_for_nonunique_equations()
        keys_for_equations_tobe_deleted = keys_for_trivial_equation_column_locations + keys_of_nonunique_equations
        # Step 4) Call equation solver "leastsq"
        N = len(self.guessvalues)
        # z1,cov,infodict,mesg,ier = leastsq(_unique_equations_to_solve, self.guessvalues,
        #                                     args=(keys_for_equations_tobe_deleted),
        #                                     maxfev= 1000*(N+2),full_output=True)
        # print("\nTotal calls made by solver:",infodict["nfev"])
        # print("result of last func call made :",infodict["fvec"])
        z = fsolve(_unique_equations_to_solve, self.guessvalues,
                   args=(keys_for_equations_tobe_deleted),
                   maxfev= 1000*(N+5))

        # CHECK SOLUTION VALIDITY
        # ------------------------
        # 1) (ΣF)in + (ΣF)out = 0
        tol = 0.01
        sumof_overallflowrates = self.npflowrates.sum(axis=0)
        # print("sumof_overallflowrates",sumof_overallflowrates)
        # print("self.flowrates", self.npflowrates)
        # 2) Create the species/component balance: (ΣFx)in - (ΣFx)out = 0
        sumof_componentflowrates = (self.npflowrates * self.npfractions).sum(axis=0)
        # print("sumof_componentflowrates",sumof_componentflowrates)
        # 3) Create the component fraction balance: (Σx) = 1
        sumof_componentfractions_ineachstream = 1-(self.npfractions.sum(axis=1))
        # print("fraction balance =", sumof_componentfractions_ineachstream)
        # 4) Perform balance-checks
        if ( np.all(np.absolute(sumof_overallflowrates) < tol) and
              np.all(np.absolute(sumof_componentflowrates) < tol) and
              np.all(np.absolute(sumof_componentfractions_ineachstream) < tol)):
            print("Unknowns successfully computed:")
            print("-------------------------------")
            print("Streams with new data have been created.")
            print("You can view the solved system by printing the returned object.\n")
            sol = PhysicalProcess("sol")._create_solvedsystem_summary( self.processname, self.npstreamnames.tolist(), self.npflowrates.tolist(),
                                                self.npfractions.tolist(), self.extrainfos)
        else:
            print("\nUnknowns could not be successfully computed")
            print("-------------------------------------------")
            print("Please check that you have correctly set up the system.")
            print("Also check your system with 'degree of freedom' analysis")
            sol = None
        return sol
   
if __name__ == "__main__":
    Mixer_NaOH = PhysicalProcess("Mixer_NaOH")
    Mixer_NaOH.Water = Mixer_NaOH.attachstreams(streamnames=["Water"])
    Mixer_NaOH.Water.setflow("F")
    Mixer_NaOH.Water.setfractions([1,0])
    Mixer_NaOH.Water.setextrainfo(["Water=0.9*Product"])
    Mixer_NaOH.NaOH = Mixer_NaOH.attachstreams(streamnames=["NaOH"])
    Mixer_NaOH.NaOH.setflow(1000)
    Mixer_NaOH.NaOH.setfractions([0,1])
    Mixer_NaOH.Product = Mixer_NaOH.attachstreams(streamnames=["Product"])
    Mixer_NaOH.Product.setflow("-F")
    Mixer_NaOH.Product.setfractions([0.9, "x"])
    print(Mixer_NaOH)
    sol2 = Mixer_NaOH.solvesystem()
    print(sol2)
    
    
    ### Extraction Distillation
    ###########################
    # https://www.youtube.com/watch?v=my1ZTIDSMbs
    # Define processes
    extraction = PhysicalProcess("Liquid-Liquid Extraction")
    distillation = PhysicalProcess("Distillation")
    overall = PhysicalProcess("Overall") # dashed boundary (a) in figure
    
    # Assign streams to processes
    extraction.m1_in = extraction.attachstreams(streamnames=["m1_in"])
    extraction.m1_in.setflow("+F")
    # Component-1: Acetic Acid, component-2: Water, component-3: Hexanol
    extraction.m1_in.setfractions([0.18, 0.82, 0])
    
    extraction.m2_in = extraction.attachstreams(streamnames=["m2_in"])
    extraction.m2_in.setflow("+F")
    extraction.m2_in.setfractions([0, 0, 1.00])
    
    extraction.m3_out = extraction.attachstreams(streamnames=["m3_out"])
    extraction.m3_out.setflow("-F")
    extraction.m3_out.setfractions(["x", 0, "x"])
    
    extraction.m4_out = extraction.attachstreams(streamnames=["m4_out"])
    extraction.m4_out.setflow("-F")
    extraction.m4_out.setfractions([0.005, 0.995, 0])
    
    distillation.m3_in = distillation.attachstreams(streamnames=["m3_in"])
    distillation.m3_in.setequalto(extraction.m3_out, flowdirection="+")
    
    distillation.m5_out = distillation.attachstreams(streamnames=["m5_out"])
    distillation.m5_out.setflow("-F")
    distillation.m5_out.setfractions([0.96, 0, 0.04])
    
    distillation.m6_out = distillation.attachstreams(streamnames=["m6_out"])
    distillation.m6_out.setflow("-F")
    distillation.m6_out.setfractions([0.028, 0, 0.972])
    distillation.m6_out.setextrainfo(["3:m6_out=0.95*m2_in"])
    
    overall.m1_in = overall.attachstreams(streamnames=["m1_in"])
    overall.m2_in = overall.attachstreams(streamnames=["m2_in"])   
    overall.m4_out = overall.attachstreams(streamnames=["m4_out"])   
    overall.m5_out = overall.attachstreams(streamnames=["m5_out"])   
    overall.m6_out = overall.attachstreams(streamnames=["m6_out"])                      
    overall.m1_in.setequalto(extraction.m1_in, flowdirection="+")
    overall.m2_in.setequalto(extraction.m2_in, flowdirection="+")
    overall.m4_out.setequalto(extraction.m4_out, flowdirection="-")
    overall.m5_out.setequalto(distillation.m5_out, flowdirection="-")
    overall.m6_out.setequalto(distillation.m6_out, flowdirection="-")
    
    # Print streams to verify they are defined correctly
    print(extraction)
    print(distillation)
    print(overall)
    
    # Perform degree of freedom analysis
    extraction.degreesoffreedom()
    distillation.degreesoffreedom()
    overall.degreesoffreedom()
    
    # From degree of freedom analysis it is clear that none
    # of the processes have degree of freedom = 0.
    # Therefore, as such the system cannot be solved.
    # However, for the 'overall; process the degree of freedom = 1.
    # If one of the flow rate is assumed, this will allow the
    # process to be solved.
    # Thus let m2 = 100
    
    overall.m2_in.setflow(100)
    
    # Check degrees of freedom again
    overall.degreesoffreedom()
    
    # The analysis now shows that "overall" system can be solved.
    overall_solution = overall.solvesystem()
    print(overall_solution)
    
    # All streams except m3 are now known.
    # To find stream m3, update information for streams
    # computed from "overall" balance.
    # Update the extraction streams
    extraction.m1_in.setequalto(overall_solution.m1_in, flowdirection="+")
    extraction.m2_in.setequalto(overall_solution.m2_in, flowdirection="+")
    extraction.m4_out.setequalto(overall_solution.m4_out, flowdirection="-")
    # Check system streams
    print(extraction)
    # Check degrees of freedom again
    extraction.degreesoffreedom()
    # a = extraction.find_unknownflowrate()
    # extraction.update_streamflowrates(a)
    sol = extraction.solvesystem()
    print(sol)

    ##### Unit1 Unit 2
    ## https://www.youtube.com/watch?v=JfD5iyoKD8w
    overall = PhysicalProcess("overall")
    mixer = PhysicalProcess("mixer")
    unit1 = PhysicalProcess("unit1")
    unit2 = PhysicalProcess("unit2")
    # attach streams to mixer
    streamnames = ["feed1", "feed2", "feedmix"]
    flowrates = [100, 10, "-F"]
    fractions = [[0.5, 0.5], [0.4, 0.6], ["x", "x"]]
    mixer.feed1, mixer.feed2, mixer.feedmix = mixer.attachstreams(streamnames=streamnames, flowrates=flowrates, fractions=fractions)
    # attach streams to unit 1
    streamnames = ["unit1feed", "top1", "bottom1"]
    flowrates = ["F", "-F", -70]
    fractions = [ ["x", "x"], [0.7, 0.3], ["x", "x"]]
    unit1.unit1feed, unit1.top1, unit1.bottom1 = unit1.attachstreams(streamnames=streamnames, flowrates=flowrates, fractions=fractions)
    # attach streams to unit 2
    streamnames = ["unit2feed", "top2", "bottom2"]
    flowrates = [70, "-F", "-F"]
    fractions = [ ["x", "x"], [0.5, 0.5], ["x", "x"]]
    unit2.unit2feed, unit2.top2, unit2.bottom2 = unit2.attachstreams(streamnames=streamnames, flowrates=flowrates, fractions=fractions)
    # attach streams to overall    
    streamnames = ["feed1", "feed2", "top1", "top2", "bottom2"]
    overall.feed1, overall.feed2, overall.top1, overall.top2, overall.bottom2 = overall.attachstreams(streamnames=streamnames)
    overall.feed1.setequalto(mixer.feed1, flowdirection="+")
    overall.feed2.setequalto(mixer.feed2, flowdirection="+")
    overall.top1.setequalto(unit1.top1, flowdirection="-")
    overall.top2.setequalto(unit2.top2, flowdirection="-")
    overall.bottom2.setequalto(unit2.bottom2, flowdirection="-")
    overall.feed2.setextrainfo(["top1=1*top2"])
   
    # print the systems to verify
    print(mixer)
    print(unit1)
    print(unit2)
    print(overall)
    
    # Compute degrees of freedom
    mixer.degreesoffreedom()
    unit1.degreesoffreedom()
    unit2.degreesoffreedom()
    overall.degreesoffreedom()
    
    # Based on Degrees of freedom analysis
    # mixer can be solved
    mixersol = mixer.solvesystem()
    print(mixersol)
    
    # Assign mixer outlet to unit1
    unit1.unit1feed.setequalto(mixersol.feedmix, flowdirection="+")
    print(unit1)
    
    
    #Find degrees of freedom again
    unit1.degreesoffreedom()
    unit2.degreesoffreedom()
    overall.degreesoffreedom()
    # DOF for unit1 = 0 ---> so solve unit1
    unit1sol = unit1.solvesystem()
    print(unit1sol)
    # Assign newly found streams
    overall.top1.setequalto(unit1sol.top1, flowdirection="-")
    unit2.unit2feed.setequalto(unit1sol.bottom1, flowdirection="+")
    print (unit2)
    print(overall)
    
    #Find degrees of freedom again
    unit2.degreesoffreedom()
    overall.degreesoffreedom()
    # DOF of overall = 0 ---> so solve overall
    overallsol = overall.solvesystem()
    print(overallsol)

    
  