"""
Created on Sat Apr 24 13:45:09 2021

@author: Harvinder Singh Gill
"""

import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import leastsq
import logging
import copy
from pychemengg.utils import multipledispatch as md

__all__ = ["PhysicalProcess"]
#---------------------------------
# Main outer class 
#---------------------------------

class PhysicalProcess(object):
    """
    A class to represent a physical process such as a mixer, 
    distillation column etc.
    
    This physical process will have streams attached to it.
    A mass balance will be performed on the streams.
    
    Usage examples: 
        mixer = PhysicalProcess("Mixer1")
        distColn1 = PhysicalProcess("Distillation Column 1")
    """  
    def __init__(self, processname):
        self.streamregister={}
        self.processname = processname

        
    #---------------------------------------
    # Inner class and its constructor method
    #---------------------------------------  
    #---------------------------------------
    # 'Stream' class constructor method
    #---------------------------------------  
    @md.multiple_dispatch(dispatch_on="kwarg_type", option="all")
    def createstreams(self, *, streamname="a_name"):
        """ version #1
        This method is called to create streams.
        Args:
            self: instance of PhysicalProcess
            streamname: name of stream to create
            Return: the stream object is returned to caller (user who is using the module)
        """           
        return PhysicalProcess.Stream(self, streamname) # return stream object to user of module
    
    @md.multiple_dispatch(dispatch_on="kwarg_type", option="all")
    def createstreams(self, *, streamnames=[[]], flowrrates=[[]], fractions=[[]], extrainfos=[[]]):
        """ version #2
        This method is called to create streams, their flowrrates
                                        their fractions, and their extrainfo.
        Args:
            self: instance of PhysicalProcess
            streamname: name of stream to create        
        """
        # if extrainfos are None for all streams, then
        # last argument when calling the function 'createstreams' will be empty list = ---> []
        # we are making it easier for user and the ymust only use [] to signify all streams 
        # are with no extrainfo
        # but a zip function is being used below, so this list must
        # have as many [] as the number of streams
        # so we do the following before calling zip
        if len(extrainfos) == 0:
            extrainfos = [[] for i in range(len(streamnames))]
            # print(f"length of extrainfos = {len(extrainfos)}")

        stream_instances_created=[]
        for stream, flow, fracts, einfo in zip(streamnames,flowrrates,fractions,extrainfos):
            tempstream = PhysicalProcess.Stream(self, stream)
            tempstream.setflow(flow)
            tempstream.setfractions(fracts)
            if len(einfo) != 0: # if no extrainfo for stream, then do not create this attribute for the stream
                tempstream.setextrainfo(einfo)
            stream_instances_created.append(tempstream)
        return stream_instances_created 
        # The stream object are returned to caller (user who is using the module). The user may capture them if they desire
        # or can simply create streams andnot bother to capture them in variables because the streams are anyway registered
        # in the streamregister and are available to code internally



            
    #---------------------------------------
    # Stream class ... Inner class
    #---------------------------------------  
    class Stream(object):
        def __init__(self, instanceof_physical_process, streamname): # here by Python rules, self is instance of class 'Stream'
            """ 
            ARGUMENTS:
                1) instanceof_physical_process :
                    TYPE: instance of class PhysicalProcess
                    DESCRIPTION: This is the PhysicalProcess to which the
                                 created stream will be bound
                2) streamname :
                     TYPE: string
                     DESCRIPTION: This is the name of the stream to be created

            RETURNS:
                Nothing
                
                Following "attributes" are created:
                    .streamname
                    .connectedprocess --> stores process linked to the stream
                    .connectedprocess.streamregister[streamname] --> stores 'Stream' object
                    
                    
        
            SIDE EFFECTS:

                
            EXCEPTIONS RAISED:

                
            RESTRICTIONS ON WHEN IT CAN BE CALLED:

            
            USAGE EXAMPLE:
 
            """
            self.streamname = streamname
            self.connectedprocess = instanceof_physical_process
            self.register_stream_with_connectedprocess(self.streamname) 

        
        # Method of inner 'Stream' class
        #---------------------------------
        def register_stream_with_connectedprocess(self,streamname): # register stream in 'streams' dictionary of the physical process
            self.connectedprocess.streamregister[streamname] = self  # here self is instance of stream passed as arg to be stored

        # Method of inner 'Stream' class
        #---------------------------------
        def setflow(self, flow): # setflow for a stream
            self.flowrate = flow            

        # Method of inner 'Stream' class
        #---------------------------------
        def setfractions(self, fracts): # set fractions of a stream
            self.fractions = fracts
            
        # Method of inner 'Stream' class
        #---------------------------------
        def setextrainfo(self, info): # set extra info of a stream
            self.extrainfo = info
        
        # Method of inner 'Stream' class
        #---------------------------------
        def deleteextrainfo(self): # delete extra info of a stream
            if hasattr(self, "extrainfo"):
                delattr(self, "extrainfo") 
        
        # Method of inner 'Stream' class
        #---------------------------------
        # copy one stream to another
        def setequalto(self, streamtocopy,  flowdirection="+"):
    
            """
            sets first argument 'self' equal to second argument
            
            """
            tempstream = copy.deepcopy(streamtocopy) # this is done so that we do not change original
                                                     # argument passed to the method
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
        def __pos__(self): # unary operator '+' performed on flow rate of a passed stream
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
            """ Method to add two streams
            1) flowrrates get added by first applying 'abs' on flowrrates
            2) component fractions are computed for the mixed stream
            """
            tempprocess = PhysicalProcess("tempprocess") # create instance of class
            tempstream = tempprocess.createstreams("temp") # create a temporary stream to return
            try: 
                totalflow = abs(self.flowrate) + abs(other.flowrate)
            except Exception:
                print("\nUnable to add stream flowrrates.")
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
                    print("Unknown fractions or flowrrates may exist in one of the streams.\n")
                    logging.exception("Caught an error in method 'add' of Stream class while computing tempfractions")
            if len(tempfractions) > 0:
                print("fractions", tempfractions)        
            tempstream.setfractions(tempfractions) # set fractions of temp stream
            del tempprocess # delete "tempprocess" 
            return tempstream # return temp stream 

        # Method of inner 'Stream' class
        #---------------------------------        
        def __sub__(self, other): # overload binary '-' operator on object 'Stream'
            """ Method to add two streams
            1) flowrrates get subtracted by first applying 'abs' on flowrrates
            2) component fractions are computed for the mixed stream
            """
            tempprocess = PhysicalProcess("tempprocess") # create instance of class to s
            tempstream = tempprocess.createstreams("temp") # create a temporary stream to return
            try: 
                totalflow = abs(self.flowrate) - abs(other.flowrate)
            except Exception:
                print("\nUnable to add stream flowrrates.")
                print("Please check flow rate values to ensure they are all known.")
                logging.exception("Caught an error in method 'add' of Stream class while computing totalflow")
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
                    print("Unknown fractions or flowrrates may exist in one of the streams.\n")
                    logging.exception("Caught an error in method 'add' of Stream class while computing tempfractions")
            if len(tempfractions) > 0:
                print("fractions", tempfractions)        
            tempstream.setfractions(tempfractions) # set fractions of temp stream
            del tempprocess # delete "tempprocess" 
            return tempstream # return temp stream 

        # Method of inner 'Stream' class
        #---------------------------------         
        def _calc_missing_xfraction(self):
            """
            Computes the missing 'x' fraction in a 
            given stream by using Σx = 1

            RETURNS
            -------
            missing 'x' == float ; or None.

            """
            if self.fractions.count("x") == 1.0:
                missingx = (1 - sum([val for val in self.fractions if val != "x"])) # because Σx = 1
                return [True, missingx]
            if self.fractions.count("x") > 1.0:
                text = "More than one 'x' unknown fractions in the stream. Try the 'solvesystem' method."
                return [False, text]
            if self.fractions.count("x") == 0.0:
                text = "There are no missing fractions."
                return [False, text]   

        # Method of inner 'Stream' class
        #---------------------------------         
        def __str__(self):

                if hasattr(self,"extrainfo"):
                    storedextrainfo = self.extrainfo
                else:
                    storedextrainfo = "None"
                return('Stream: {streamname}\n'
                'Attached to physical process: {process.processname}\n'
                'flowrate: {flowrate}\n'
                'fractions: {fractions}\n'
                ).format(**self.__dict__) + "extrainfo: " + storedextrainfo
           
    #-----------------------------------------
    # Methods for outer class "PhysicalProcess"
    #-----------------------------------------
    
    # Method for outer class "PhysicalProcess"
    #----------------------------------------    
    def _make_numpyarrayof_fractions(self):
        a = [streamobject.fractions for streamname, streamobject in self.streamregister.items()]
        return np.array(a)
    
    # Method for outer class "PhysicalProcess"
    #----------------------------------------                
    def _make_numpyarrayof_flowrates(self):
        a = [streamobject.flowrate for streamname, streamobject in self.streamregister.items()]
        # print("makenumpyFlowrate = ", np.vstack(a))
        return np.vstack(a)

    # Method for outer class "PhysicalProcess"
    #----------------------------------------    
    def _make_numpyarrayof_streamnames(self):
        a = [streamname for streamname, streamobject in self.streamregister.items()]
        return np.vstack(a)
    
    # Method for outer class "PhysicalProcess"
    #----------------------------------------    
    def _make_arrayof_extrainfo(self):
        a=[]
        for streamname, streamobject in self.streamregister.items():
            if hasattr(streamobject, 'extrainfo'):
                a.append(streamobject.extrainfo)
            else:
                a.append([])
        return a
    
    # Method for outer class "PhysicalProcess"
    #-----------------------------------------    
    def _parse_extrainfo(self):
        """
        This method parses the extra info for 'TOTAL flowrate' relationship OR
        'COMPONENT flowrate' relationships between streams.

        Returns
        -------
        parsed_extrainfo : as list of lists
                            [[], [],[]]
            DESCRIPTION: Each sub list contains the following list elements 
          
            CASE 1: COMPONENT flowrate relationship
            [componentindex, left_streamsign, left_streamname, ratio, right_streamname]

            example:
                let extrainfo : ["2:S1=1.2*S4"] 
                This means for component '2', the component flowrate in S1 is 1.2 times that in S4
                So component '2' flowrrates in S1 and S4 will be
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
        for streamname, streamobject in self.streamregister.items():
            if hasattr(streamobject, 'extrainfo'):
                for items in streamobject.extrainfo:
                    # first parse for 'component flowrate'
                    # if info is not for 'component flowrate'
                    # the first statement of 'try' will fail
                    # it then helps to decide that the relationship
                    # is for total stream flowrrates
                    try:
                        componentindex, other_1 = items.split(":")
                        left_streamname, other_2 = other_1.split("=")
                        ratio, right_streamname = other_2.split("*")
                        left_streamsign = self._getsignof(self.streamregister[left_streamname].flowrate)
                        right_streamsign = self._getsignof(self.streamregister[right_streamname].flowrate)
                        temp = [componentindex, left_streamsign, left_streamname, ratio, right_streamsign, right_streamname]
                        parsed_extrainfo.append(temp)
                    except ValueError:
                        # print(items)
                        left_streamname, other_2 = items.split("=")
                        ratio, right_streamname = other_2.split("*")
                        left_streamsign = self._getsignof(self.streamregister[left_streamname].flowrate)
                        right_streamsign = self._getsignof(self.streamregister[right_streamname].flowrate)
                        temp = [False, left_streamsign, left_streamname, ratio, right_streamsign, right_streamname]
                        # print("extra Info", temp)
                        parsed_extrainfo.append(temp)
        return parsed_extrainfo


    # Method for outer class "PhysicalProcess"
    #-----------------------------------------
    
    def _make_numpyarraysof_streamdata(self):
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
        # arg : Input parameter = list of flowrrates
        a = np.where(arg=="F")
        return a
    
    # Method for outer class "PhysicalProcess"
    #----------------------------------------      
    @staticmethod
    def _find_whereminusFs(arg):
        # arg : Input parameter = list of flowrrates
        a = np.where(arg=="-F")
        return a
    
    # Method for outer class "PhysicalProcess"
    #----------------------------------------      
    @staticmethod
    def _find_wherexs(arg):
        # arg : Input parameter = list of fractions
        #                         of single stream
        a = np.where(arg=="x")
        return a
   
    # Method for outer class "PhysicalProcess"
    #-----------------------------------------   
    def __str__(self, print_thesestreams=None):
        if print_thesestreams is None:
            streamstoprint = self.streamregister
        if print_thesestreams is not None:
            streamstoprint = print_thesestreams            
        # print header info before printing data
        headerstring = "Process streams for '" + self.processname + "' are :"
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
        maxflow = max([len(str(round(vals.flowrate,2))) for vals in streamstoprint.values() if type(vals.flowrate) != str], default=0)
        
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
        # Case: Apply to all streams
        unknownx_detailsfound = []
        if streamname == None:
            for strname, streamobject in self.streamregister.items():
                missingx = streamobject._calc_missing_xfraction()
                if missingx[0] == True:
                    location_wherexunknown = [i for i, val in enumerate(streamobject.fractions) if val =="x"]
                    unknownx_detailsfound.append([True, strname, location_wherexunknown[0], missingx[1]])
                if missingx[0] == False:
                    unknownx_detailsfound.append([False, strname, missingx[1]])

        # Case: Apply balance to single stream when 'streamname' provided in the method argument
        elif streamname != None:
            streamobject = self.streamregister[streamname]
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
        """
        This method will update stream fraction at ONE given component-location
        
        new_fractioninfo argument has the form:
        
            [True, streamname, location_wherexunknown, missingx]
        or [False, streamname, missingx]
        
        """
        for info in new_fractioninfo:
            if info[0] == True:
                self.streamregister[info[1]].fractions[info[2]]=info[3]
                print("For ", '"', self.processname, '"', ": the stream", '"', info[1], '"', "The new component fraction @ position =", '"',info[2],'"', "is now = ", info[3])
                
            if info[0] == False:
                print("For ", '"', self.processname, '"', ": the stream", '"', info[1], '"', info[2])
 
    # Method for outer class "PhysicalProcess"
    #----------------------------------------- 
    def find_unknownflowrate(self): # ΣFin + ΣFout = 0
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
            text = "More than one 'F'/'-F' unknown flowrrates in stream. Try the 'solvesystem' method"
            unknownflowrate_detailsfound.append([False, text])
            return unknownflowrate_detailsfound
        
        if (unknownplusF_count + unknownminusF_count) == 0.0:
            text = "There are no missing flowrrates"
            unknownflowrate_detailsfound.append([False, text])
            return unknownflowrate_detailsfound


    # Method for outer class "PhysicalProcess"
    #--------------------------------------------     
    def update_streamflowrates(self, newflow):
        # print(newflow)
        for streamdata in newflow:
            if streamdata[0] == True:
                self.streamregister[streamdata[1]].flowrate=streamdata[2]
                print("For ", '"', self.processname, '"', ": the stream",'"', streamdata[1], '"', "has new flowrate =", streamdata[2])

            if streamdata[0] == False:
                print("For ", '"', self.processname, '"', streamdata[1])
                


    # Method for outer class "PhysicalProcess"
    #-----------------------------------------
    def perform_component_massbalance(self):
        # Performs (ΣFx )in = (ΣFx)out for each component
        # Returns numpy array of mass balance for each component
        # if the process is uniquely defined then these terms will
        # be close to zero
        self._make_numpyarraysof_streamdata()
        print(self.npflowrates, self.npfractions)
        npcomponent_massbalance = (self.npflowrates * self.npfractions).sum(axis=0)
        
        # THIS IS ALTERNATE CODE
        #------------------------
        # npstream_flowandfraction_data = [[np.array(streamobject.flowrate),np.array(streamobject.fractions)]
        #                for streamname, streamobject in self.streamregister.items()]
        # component_mass = [numpyFlow * numpyFractions for numpyFlow, numpyFractions in npstream_flowandfraction_data]
        # npcomponent_mass = np.array(component_mass)
        # npcomponent_massbalance = npcomponent_mass.sum(axis=0) 
        
        return npcomponent_massbalance # return type <class 'numpy.ndarray'> ... 1 x n array
    
    
    # Method for outer class "PhysicalProcess"
    #----------------------------------------- 
    def _find_unknownxandFlocations(self):
        self._make_numpyarraysof_streamdata()
        plusF_location = self._find_whereplusFs(self.npflowrates)
        minusF_location = self._find_whereminusFs(self.npflowrates)
        x_location = self._find_wherexs(self.npfractions)
        return (plusF_location, minusF_location, x_location)  
    
    
    # Method for outer class "PhysicalProcess"
    #-----------------------------------------   
    def _getsignof(self, F):
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
        self.totalequations_possible = self.extrainfo_count + self.component_count + self.streamswith_unknownx_count
   
        print("Degrees of freedom analysis")
        print("===========================")
        print("Number of unknown flowrrates :-->", self.unknown_Fs_count)
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
            return None
        
        if self.totalunknowns == self.totalequations_possible:
            print("\nSystem can be solved")
            print("\nThere are = " , self.totalunknowns , "unknowns and" , self.totalequations_possible , "possible equations")
            return None

    # Method for outer class "PhysicalProcess"
    #-----------------------------------------           
    def solvesystem(self):
        
        self.degreesoffreedom() # perform degrees of freedom analysis
       
        # Inner method for outer class "PhysicalProcess"
        #-----------------------------------------------
        def _solverfunction(guess):
            # This function is called by fsolve/leastsq.
            # It returns the equations of the physical process as F[]
            # Input parameter = Initial guess value
            # print("guess = ", guess)

            
            if self.totalunknowns > self.totalequations_possible:
                print("\nSystem is underspecified because there are more unknowns than available equations")
                print("\nThere are = ", self.totalunknowns, "unknowns and", self.totalequations_possible,
                      "equations")
                return None
            
            if self.totalunknowns <= self.totalequations_possible:
                guessindex = 0 # initialize index to array of guess value
                
                # (A) ASSIGN GUESS VALUES
                # -----------------------
                # 1) First assign 'guess' values to unknown 'x' 
                if self.unknownx_count > 0:
                    self.npfractions[(self.x_location)] = guess[guessindex : self.unknownx_count]
                    guessindex = guessindex + self.unknownx_count
                self.npfractions = self.npfractions.astype(float)
                # 2) Assign 'guess' values to unknown 'F'
                if self.unknownplusF_count > 0:
                    self.npflowrates[(self.plusF_location)] = guess[guessindex : guessindex + self.unknownplusF_count]
                    guessindex = guessindex + self.unknownplusF_count
                # 3) Assign 'guess' values to unknown '-F'
                if self.unknownminusF_count > 0:
                    self.npflowrates[(self.minusF_location)] = guess[guessindex : guessindex + self.unknownminusF_count]
                    guessindex = guessindex + self.unknownminusF_count
                self.npflowrates = self.npflowrates.astype(float)
                
                # (B) NEXT CREATE EQUATIONS TO ALLOW DETERMINATION OF 'x' and 'F'/'-F'
                # --------------------------------------------------------------------
                # 1) Start by forming equations using'Extra Info' of streams
                if self.extrainfo_count > 0:
                    extrainfoindex = 0 # initialize index into extrainfo relations
                    extrarelations ={}
                    for indx, relation in enumerate(self.extrainfo_records):
                        isextrainfo_between_massofstreams = relation[0]
                        left_streamsign = relation[1]
                        left_streamname = relation[2]
                        ratio = relation[3]
                        right_streamsign = relation[4]
                        right_streamname = relation[5]
                        left_streamflow = self.npflowrates[np.where(self.npstreamnames == left_streamname)]
                        right_streamflow = self.npflowrates[np.where(self.npstreamnames == right_streamname)]
                        if isextrainfo_between_massofstreams == False:
                            # print("leftflow:", left_streamflow, " rightflow: ", right_streamflow)
                            if left_streamsign * right_streamsign > 0: # both streams are either '+':inlet or '-':outlet
                                extrarelations[extrainfoindex] = left_streamsign*abs(left_streamflow) - right_streamsign * abs(((eval(ratio) * right_streamflow)))
                                # print(">0", left_streamflow - ((eval(ratio) * right_streamflow)))
                                extrainfoindex = extrainfoindex + 1
                            if left_streamsign * right_streamsign < 0: # one stream is '+':inlet other is '-':outlet
                                extrarelations[extrainfoindex] = (left_streamsign*(left_streamflow)) - abs(right_streamsign * (((eval(ratio) * right_streamflow))))
                                # extrarelations[extrainfoindex] = left_streamsign*abs(left_streamflow) + right_streamsign * abs(((eval(ratio) * right_streamflow)))
                                # print("<0", left_streamflow + ((eval(ratio) * right_streamflow)))
                                extrainfoindex = extrainfoindex + 1
                        if isextrainfo_between_massofstreams != False: # then extra relationship is between total stream flowrrates
                            column_coordinate = int(isextrainfo_between_massofstreams[0])-1
                            left_rowcoordinate = (np.where(self.npstreamnames == left_streamname))
                            right_rowcoordinate = (np.where(self.npstreamnames == right_streamname)) 
                            left_streamfraction = self.npfractions[left_rowcoordinate[0][0], column_coordinate]
                            right_streamfraction = self.npfractions[right_rowcoordinate[0][0], column_coordinate]
                            if left_streamsign * right_streamsign > 0: # both streams are either '+':inlet or '-':outlet
                                extrarelations[extrainfoindex] = left_streamsign*abs(left_streamfraction * left_streamflow) - right_streamsign * abs(eval(ratio) * right_streamfraction * right_streamflow) 
                                extrainfoindex = extrainfoindex + 1
                            if left_streamsign * right_streamsign < 0: # one stream is '+':inlet other is '-':outlet
                                extrarelations[extrainfoindex] = left_streamsign*(left_streamfraction * left_streamflow) - abs(right_streamsign * eval(ratio) * right_streamfraction * right_streamflow) 
                                extrainfoindex = extrainfoindex + 1
                # 2) Create the component/species balance: (ΣFx)in - (ΣFx)out = 0
                componentmassbalance = (self.npflowrates * self.npfractions).sum(axis=0)
                # 3) Create the component fraction balance: (Σx) = 0
                componentmassfractionbalance = 1-(self.npfractions.sum(axis=1))
                
                # (C) NEXT DEFINE FUNCTION THAT IS ASSIGNED DIFFERENT EQUATIONS FROM ABOVE
                #     THIS FUNCTION IS THEN RETURNED TO 'fsolve'
                # ------------------------------------------------------------------------
                var_equations = np.ones(len(guess)) # These are equations that will contain 'x' and '+F'/'-F' and 
                                                    # will be solved to obtain unknown quantities
                vareqindex = 0 # 
                # 1) Start by assigning 'EXTRA INFO' equations to var_equations        
                if self.extrainfo_count > 0:
                    for equation in extrarelations.values():
                        # print("Inside extra info ... var_equations-> allocation")
                        var_equations[vareqindex] = equation
                        vareqindex = vareqindex + 1
                        # print("vareqindex = ", vareqindex)
                        if vareqindex >= self.totalunknowns:
                            break
                # 2) Check if enough var_equations[] have been formed
                # If not, then use SPECIES balances to form var_equations
                # but avoid using the balances if in that equation
                # all 'Fs' and 'xs' are known                
                if ( vareqindex ) < self.totalunknowns: 
                    if self.unknownx_count > 0:
                        arraycolumns_withxunknown = set( self.x_location[1] )
                    elif self.unknownx_count <=0:
                        arraycolumns_withxunknown = set(  )
                    for column, val in enumerate(componentmassbalance):
                        # the 'or' condition below means to take component
                        # balance equation even if all 'x' are known for that
                        # component for all streams, but at least one 'F'/'-F'
                        # is not known. This is because ΣFixi = 0
                        if ( column in arraycolumns_withxunknown ) or ( self.unknown_Fs_count > 0 ):
                            var_equations[vareqindex] = val
                            vareqindex = vareqindex + 1
                            # print("inside component balance")
                            # print("x ", val)
                            if (vareqindex) >= self.totalunknowns:
                                break
                # NOTE: 
                # The following 'OVERALL FLOW BALANCE' was REMOVED FROM CODE LOGIC
                # ----------------------------------------------------------------
                # check if enough var_equations[] have been formed     
                # Use OVERALL FLOW BALANCE as an equation if 'F' or '-F' is an unknown
                # Overall balance is of no use if all flowrrates are known
                # Furthermore we select overall flow balance if more than 
                # one flowrate is unknown
                # NOTE: Since Overall Balance is derivable from Species balance
                # when this equation was selected, in certain cases the solutions could
                # not be found. This is because the Overall and Species balances are dependent
                # So this was taken out from the logic :
                # if ( vareqindex ) < totalunknowns:
                #     if ( unknown_Fs_count > 1 ):
                #         # print("Inside OVERALL BALANCE ... info F-> allocation")
                #         F[vareqindex] = self.npflowrates.sum(axis = 0)
                #         vareqindex = vareqindex + 1
                #         # print("vareqindex = ", vareqindex)
                         
                
                # 3) Check if enough var_equations[] have been formed
                # if not, Use FRACTION balances
                # but don't use streams for which all 'x' are known
                if (vareqindex) < self.totalunknowns:
                    arrayrows_withxunknown = set(self.x_location[0])
                    for row, val in enumerate(componentmassfractionbalance):
                        if row in arrayrows_withxunknown:
                            var_equations[vareqindex] = val
                            vareqindex = vareqindex + 1
                            # print("in fraction balance")
                            if (vareqindex) >= self.totalunknowns:
                                break   
            return var_equations
        
        # PREPARE GUESS VALUES TO CALL 'fsolve'
        guessvalues = np.ones(self.totalunknowns)
        guessindex = 0
        # guessvalues are used to assign values to unknown 'x' and 'F'
        # First the 'x' values are assigned, and then 'F' values
        # So we construct guessvalues such that it's initial elements
        # corresponding to 'x' are initialized to 0.5 
        # since 'x' lies between 0 and 1
        if self.unknown_x_count > 0:
            for i in range(self.unknown_x_count):
                guessvalues[guessindex] = 0.5
                guessindex = guessindex + 1
        
        # The remainder elements of guessvalues must be set equal to flow rates
        # that are already known
        # So we find known flow rate values to use them in initilization
        # If we randomly initialize guess values of flowrrates to say 1, 100, or 1000 etc
        # then the fsolve function maynot converge because these random values
        # maybe far from the true solution
        knownflowrates = [float(flow[0]) for flow in self.npflowrates if (flow[0] != "F" and flow[0] != "-F")]
        sorted_knownflowrates = sorted(knownflowrates, reverse=True)
        if self.unknown_Fs_count > 0:
            for i in range(self.unknown_Fs_count):
                if i < len(knownflowrates):
                    guessvalues[guessindex] = sorted_knownflowrates[i]
                    guessindex = guessindex + 1
                else: # use the last flowrate in list as guess value
                    guessvalues[guessindex] = knownflowrates[-1]
                    guessindex = guessindex + 1
        
        # Call the '_solverfunction' to solve and find solution
        z = fsolve( _solverfunction, guessvalues )
        
        # Check Solution
        # --------------
        # 1) (ΣF)in + (ΣFx)out = 0
        tol = 0.01
        sumof_overallflowrates = self.npflowrates.sum(axis=0)
        # 2) Creates the species/component balance: (ΣFx)in - (ΣFx)out = 0
        sumof_componentflowrates = (self.npflowrates * self.npfractions).sum(axis=0)
        # print("component balance = ", sumof_componentflowrates)
        # 3) Create the component fraction balance: (Σx) = 0
        sumof_componentfractions_ineachstream = 1-(self.npfractions.sum(axis=1))
        # print("fraction balance =", sumof_componentfractions_ineachstream)
        # 4) Perform balance checks
        if ( np.all(np.absolute(sumof_overallflowrates) < tol) and
             np.all(np.absolute(sumof_componentflowrates) < tol) and
             np.all(np.absolute(sumof_componentfractions_ineachstream) < tol)):
            print("\nUnknowns successfully computed:")
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
   
        
    # Method for outer class "PhysicalProcess"
    #---------------------------------    
    @staticmethod
    def _create_solvedsystem_summary(processname, npstreamnames, npflowrates, npfractions, extrainfos ):
        string = '"Solution" to ' + processname
        solvedsystem = PhysicalProcess(string)
        for num, strmName in enumerate(npstreamnames):
            # headerstring = "solvedsystem"+"."+strmName[0]+"="+"solvedsystem.createstreams"+"("+'"'+strmName[0]+'"'+" )"
            headerstring = "solvedsystem."+strmName[0]+"= solvedsystem.createstreams(streamname="+"""strmName[0]"""+" )"
            # headerstring = "solvedsystem."+strmName[0]+"= solvedsystem.createstreams("+"""strmName[0]"""+" )"
            exec(headerstring)
            solvedsystem.streamregister[strmName[0]].flowrate = npflowrates[num][0]
            solvedsystem.streamregister[strmName[0]].fractions = npfractions[num]
            if len(extrainfos[num]) > 0:
                solvedsystem.streamregister[strmName[0]].extrainfo = extrainfos[num]
        return solvedsystem


if __name__ == "__main__":
##########################################
##########################################
    # Example 1
    #--------------     
    Mixer1 = PhysicalProcess("Mixer1")
    Mixer1.S1 = Mixer1.createstreams(streamname="S1")
    Mixer1.S1.setflow(500)
    Mixer1.S1.setfractions([0.2,0.8])  
    Mixer1.S2 = Mixer1.createstreams(streamname="S2")
    # S2.setflow(340)
    Mixer1.S2.setflow("F")
    Mixer1.S2.setfractions([0.45,0.55])
    Mixer1.S3 = Mixer1.createstreams(streamname="S3")
    Mixer1.S3.setflow("-F")
    # S3.setflow(-210)
    Mixer1.S3.setfractions([0.7,0.3])
    Mixer1.S4 = Mixer1.createstreams(streamname="S4")
    Mixer1.S4.setflow("-F")
    # S4.setflow(-630)
    # S4.setfractions(["x", "x"])
    Mixer1.S4.setfractions([0.168253968, 0.831746032])
    Mixer1.S4.setextrainfo(["S4=3*S3"])
    print(Mixer1) 
    sol = Mixer1.solvesystem()
    print(sol)
    # out = sol.S3 + sol.S4 - sol.S2
    # print("------------add streams S3 and S4 -----------")
    # print(out.flowrate)
    # print(out.fractions)
    
    # ##########################################
    # ##########################################    
    # # # Example 2
    # # #------------
    
    # Mixer2 = PhysicalProcess("Mixer2")
    # streams = ["S1", "S2", "S3", "S4"]
    # for strm in streams:
    #     globals()[strm] = Mixer2.createstreams(strm)
    # # S1 = Mixer2.createstreams("S1")
    # # # S1.setflow(500)
    # S1.setflow("F")
    # # S1.setextrainfo(["S1=500/340*S2"])
    # S1.setfractions([0.2, 0.16, 0.64])
    # # S2 = Mixer2.createstreams("S2")
    # S2.setflow(340)
    # # # S2.setfractions([0.45, 0.32, 0.23])
    # S2.setfractions(["x",0.32, "x"])
    # # S3 = Mixer2.createstreams("S3")
    # # # S3.setflow(-210)
    # S3.setflow("-F")
    # # # S3.setfractions([0.3, 0.3, 0.4])
    # S3.setfractions([0.3, 0.3, "x"])
    # # S4 = Mixer2.createstreams("S4")
    # S4.setflow(-630)
    # # S4.setflow("-F")
    # # # S4.setfractions([0.301587, 0.199683, 0.49873])
    # S4.setfractions([0.301587, "x", "x"])
    # S4.setextrainfo(["S4=630/210*S3"])
    # print(Mixer2)
    # # Mixer2._find_unknownxandFlocations()
    # se = Mixer2.solvesystem()
    # # print(Mixer2)
    # print(se)  
    # # out = S1+S2
    # # out.flowrate
    
    # ##########################################
    # ##########################################    
    # # # # Example 3
    # # # # --------
    # # # # Mixer1 = selfs("Mixer1", streamnames = ["Water","NaOH","Product"],
    # # # #                         stream_flowrates=["F", 1000, "-F"],
    # # # #                         stream_component_fractions=[[1, 0],
    # # # #                                                   [0, 1],
    # # # #                                                   ["x", "x"]],
    # # # #                         # extra_info = None)
    # # # #                         extra_info = ["Water=0.9*Product"])   
    
    # Mixer3 = PhysicalProcess("Mixer3")
    # Mixer3.Water = Mixer3.createstreams("Water")
    # Mixer3.Water.setflow("F")
    # # # Mixer3.Water.setflow(9000)
    # Mixer3.Water.setfractions([1,0])
    # Mixer3.Water.setextrainfo(["Water=0.9*Product"])
    # Mixer3.NaOH = Mixer3.createstreams("NaOH")
    # Mixer3.NaOH.setflow(1000)
    # Mixer3.NaOH.setfractions([0,1])
    # Mixer3.Product = Mixer3.createstreams("Product")
    # Mixer3.Product.setflow("-F")
    # # # Product.setflow(-10000)
    # Mixer3.Product.setfractions(["x", 0.1])
    # # # Product.setfractions([0.9, 0.1])
    
    # print(Mixer3)
    # sol2 = Mixer3.solvesystem()
    
    
    # print(sol2)
    # aer3 = Mixer3.find_unknownflowrate()   
    
       
    # ##########################################
    # ##########################################
    # Mixer3a = PhysicalProcess("Mixer3a")
    # Mixer3a.Water = Mixer3a.createstreams("Water")
    # # Water.setflow("F")
    # Mixer3a.Water.setflow(9000)
    # Mixer3a.Water.setfractions([1,0])
    #    # Water.setextrainfo(["Water=0.9*Product"])
    # Mixer3a.NaOH = Mixer3a.createstreams("NaOH")
    # Mixer3a.NaOH.setflow(1000)
    # Mixer3a.NaOH.setfractions([0,1])
    # Mixer3a.Product = Mixer3a.createstreams("Product")
    # Mixer3a.Product.setflow("-F")
    # # # Product.setflow(-10000)
    # Mixer3a.Product.setfractions(["x", "x"])
    # # # Product.setfractions([0.9, 0.1])
    
    # sol2a = Mixer3a.solvesystem()
    # print(Mixer3a)
    
    # print(sol2a)
    # aer3a = Mixer3a.find_unknownflowrate()
    
    # ##########################################
    # ##########################################
    # Mixer3b = PhysicalProcess("Mixer3b")
    # Mixer3b.Water = Mixer3b.createstreams("Water")
    # Mixer3b.Water.setflow("F")
    # Mixer3b.Water.setflow(9000)
    # Mixer3b.Water.setfractions([1,0])
    # # Water.setextrainfo(["Water=0.9*Product"])
    # Mixer3b.NaOH = Mixer3b.createstreams("NaOH")
    # Mixer3b.NaOH.setflow(1000)
    # Mixer3b.NaOH.setfractions([0,1])
    # Mixer3b.Product = Mixer3b.createstreams("Product")
    # Mixer3b.Product.setflow("-F")
    # # # Product.setflow(-10000)
    # Mixer3b.Product.setfractions(["x", 0.1])
    # Mixer3b.Product.setextrainfo(["Product=(1/0.9)*Water"])
    # # # Mixer3b.Product.setfractions([0.9, 0.1])
    
    # print(Mixer3b)
    # sol2 = Mixer3b.solvesystem()
    
    
    # print(sol2)
    # aer3 = Mixer3b.find_unknownflowrate()
    
    # ##########################################
    # ##########################################
    # Mixer4=PhysicalProcess("Mixer4")
    # Mixer4.S1 = Mixer4.createstreams("S1" )
    # Mixer4.S1.setflow(100)
    # # # Mixer4.S1.setflow("F")
    # Mixer4.S1.setfractions( [0.5, 0.4, 0.1] )
    # Mixer4.S2 = Mixer4.createstreams( "S2" ) 
    # # # Mixer4.S2.setflow("-F")
    # Mixer4.S2.setflow(-60)
    # Mixer4.S2.setfractions( [0.8, 0.05, 0.15] )
    # Mixer4.S3 = Mixer4.createstreams( "S3" )
    # # Mixer4.S3.setflow("-F")
    # Mixer4.S3.setflow(-40)
    # Mixer4.S3.setfractions(["x", "x","x"])
    # # # Mixer4.S3.setfractions([0.05, 0.925,  0.025])
    # # print(Mixer4)
    
    # sol = Mixer4.solvesystem()
    # print(Mixer4)
    # # print(Mixer4._getsignof("F"))
    
    # print(sol)
    
    # ##########################################
    # ##########################################    
    # Mixer5 = PhysicalProcess("Mixer5")
    # Mixer5.S1 = Mixer5.createstreams("S1")
    # Mixer5.S1.setflow(500)
    # # Mixer5.S1.setflow("F")
    # # Mixer5.S1.setextrainfo(["S1=500/340*S2"])
    # Mixer5.S1.setfractions([0.2, 0.16, 0.64])
    # Mixer5.S2 = Mixer5.createstreams("S2")
    # Mixer5.S2.setflow(340)
    # # # Mixer5.S2.setfractions([0.45, 0.32, 0.23])
    # Mixer5.S2.setfractions(["x",0.32, "x"])
    # Mixer5.S3 = Mixer5.createstreams("S3")
    # Mixer5.S3.setflow(-840)
    # # Mixer5.S3.setflow("-F")
    # # Mixer5.S3.setfractions([0.30119, 0.224762, 0.474048])
    # Mixer5.S3.setfractions(["x", 0.224762, 0.474048])
    
    # sol = Mixer5.solvesystem()
    # print(Mixer5)
    
    # print(sol)
    # aer5 = Mixer5.find_unknownflowrate()
    
    # ##########################################
    # ##########################################     
    # Mixer6=PhysicalProcess("Mixer6")
    # Mixer6.S1 = Mixer6.createstreams("S1" )
    # # Mixer6.S1.setflow(100)
    # Mixer6.S1.setflow("F")
    # Mixer6.S1.setfractions( [0.5, 0.4, 0.1] )
    # Mixer6.S2 = Mixer6.createstreams( "S2" ) 
    # # # Mixer6.S2.setflow("-F")
    # Mixer6.S2.setflow(-60)
    # # Mixer6.S2.setfractions( [0.8, 0.05, 0.15] )
    # Mixer6.S2.setfractions( ["x", 0.05, "x"] )
    # Mixer6.S3 = Mixer6.createstreams( "S3" )
    # # Mixer6.S3.setflow("-F")
    # Mixer6.S3.setflow(-40)
    # Mixer6.S3.setfractions(["x", 0.925,0.025])
    # # # Mixer6.S3.setfractions([0.05, 0.925,  0.025])
    # print(Mixer6)
    
    # ans = Mixer6.solvesystem()
    # print(ans)
    
    # # print(Mixer6._getsignof("F"))
    # aer = Mixer6.find_unknownflowrate()
    
    # ##########################################
    # ##########################################    
    # Mixer7 = PhysicalProcess("Ethanol Mixer")
    # Mixer7.S1 = Mixer7.createstreams("S1")
    # Mixer7.S2 = Mixer7.createstreams("S2")
    # Mixer7.S3 = Mixer7.createstreams("S3")
    # Mixer7.S1.setflow(-10)
    # Mixer7.S1.setfractions([0.05, 1-0.05])
    # Mixer7.S2.setflow("F")
    # Mixer7.S2.setfractions([0.01, 1-0.01])
    # Mixer7.S3.setflow("F")
    # Mixer7.S3.setfractions([0.41, 1-0.41])
    # Mixer7.sol = Mixer7.solvesystem()
    # print(Mixer7.sol)
    
    # ##########################################
    # ##########################################    
    # column12=PhysicalProcess("Distillation Column")
    # column12.S1 = column12.createstreams("S1" )
    # # column12.S1.setflow(100)
    # column12.S1.setflow("F")
    # column12.S1.setfractions( [0, 0.03, 0.97] )
    # # column12.S1.setextrainfo(["3:S5=0.90*S1"])
    # column12.S2 = column12.createstreams( "S2" ) 
    # # # column12.S2.setflow("-F")
    # column12.S2.setflow(5300)
    # # column12.S2.setfractions( [0.8, 0.05, 0.15] )
    # column12.S2.setfractions( ["x", "x", 0] )
    # column12.S3 = column12.createstreams( "S3" )
    # column12.S3.setflow("-F")
    # column12.S3.setfractions([1, 0, 0])
    # column12.S3.setextrainfo(["S3=0.5*S1"])
    # column12.S4 = column12.createstreams( "S4" )
    # column12.S4.setflow(-1200)
    # column12.S4.setfractions([0.7, "x", "x"])
    # column12.S5 = column12.createstreams( "S5" )
    # column12.S5.setflow("-F")
    # column12.S5.setfractions([0, 0.6, 0.4])
    # column12.S5.setextrainfo(["3:S5=0.90*S1"])
    # print(column12)
    
    # ans = column12.solvesystem()
    # print(ans)
    
    # ##########################################
    # ##########################################
    # feedMixer = PhysicalProcess("Feed mixer")
    # feedMixer.S1 = feedMixer.createstreams("S1")
    # feedMixer.S1.setflow(5000)
    # feedMixer.S1.setfractions( [0.2, 1-0.2] )
    # feedMixer.S2 = feedMixer.createstreams("S2")
    # feedMixer.S2.setflow("F")
    # feedMixer.S2.setfractions( [0.3, 1-0.3] )
    # feedMixer.S3 = feedMixer.createstreams("S3")
    # feedMixer.S3.setflow("-F")
    # feedMixer.S3.setfractions(["x","x"])
    # print(feedMixer)
    
    # evaporator = PhysicalProcess("Evaporator")
    # evaporator.S4 = evaporator.createstreams("S4")
    # evaporator.S4.setequalto(feedMixer.S3, flowdirection = "+")
    # # evaporator.S4.setflow("F")
    # evaporator.S5 = evaporator.createstreams("S5")
    # evaporator.S5.setflow("-F")
    # evaporator.S5.setfractions([0,1])
    # evaporator.S6 = evaporator.createstreams("S6")
    # evaporator.S6.setflow("-F")
    # evaporator.S6.setfractions([0.35, 1-0.35])
    # print(evaporator)
    
    # crystallizer = PhysicalProcess("Crystallizer")
    # crystallizer.S7 = crystallizer.createstreams("S7")
    # crystallizer.S7.setequalto(evaporator.S6, flowdirection = "+")
    # # crystallizer.S7.setflow("F")
    # crystallizer.S8 = crystallizer.createstreams("S8")
    # crystallizer.S8.setflow("-F")
    # crystallizer.S8.setfractions([1,0])
    # crystallizer.S9 = crystallizer.createstreams("S9")
    # crystallizer.S9.setflow("-F")
    # crystallizer.S9.setfractions([0.3,1-0.3])
    # print(crystallizer)
    
    # splitter = PhysicalProcess("Splitter")
    # splitter.S10 = splitter.createstreams("S10")
    # splitter.S10.setequalto(crystallizer.S9, flowdirection = "+")
    # # splitter.S10.setflow("F")
    # splitter.S11 = splitter.createstreams("S11")
    # splitter.S11.setequalto(feedMixer.streamregister["S2"], flowdirection = "-")
    # splitter.S11.setflow("-F")
    # splitter.S12 = splitter.createstreams("S12")
    # splitter.S12.setflow("-F")
    # splitter.S12.setfractions([0.3, 1-0.3])
    # print(splitter)
    
    # overall = PhysicalProcess("Overall")
    # overall.S13 = overall.createstreams("S13")
    # overall.S13.setequalto(feedMixer.S1, flowdirection="+") 
    # overall.S14 = overall.createstreams("S14")
    # overall.S14.setequalto(evaporator.S5, flowdirection="-")
    # overall.S15 = overall.createstreams("S15")
    # overall.S15.setequalto(crystallizer.S8, flowdirection="-")
    # overall.S16 = overall.createstreams("S16")
    # overall.S16.setequalto(splitter.S12, flowdirection="-")
    # print(overall)
    
    
    # crystallizer = PhysicalProcess("Crystallizer")
    # crystallizer.S7 = crystallizer.createstreams("S7")
    # crystallizer.S7.setequalto(evaporator.S6, flowdirection = "+")
    # # crystallizer.S7.setflow("F")
    # crystallizer.S8 = crystallizer.createstreams("S8")
    # crystallizer.S8.setflow("-F")
    # crystallizer.S8.setfractions([1,0])
    # crystallizer.S9 = crystallizer.createstreams("S9")
    # crystallizer.S9.setflow("-F")
    # crystallizer.S9.setfractions([0.3,1-0.3])
    # print(crystallizer)
    
    # overall = PhysicalProcess("Overall")
    # overall.S13 = overall.createstreams("S13")
    # # overall.S13.setflow(5000)
    # # overall.S13.setfractions([0.2, 0.8])
    # overall.S13.setequalto(feedMixer.S1, flowdirection="+") 
    # overall.S14 = overall.createstreams("S14")
    # # overall.S14.setflow("-F")
    # # overall.S14.setfractions([0,1])
    # overall.S14.setequalto(evaporator.S5,flowdirection="-")
    # overall.S15 = overall.createstreams("S15")
    # # overall.S15.setflow("-F")
    # # overall.S15.setfractions([1,0])
    # # overall.S15.setextrainfo(["S15=25*S16"])
    
    # overall.S15.setequalto(crystallizer.S8,flowdirection="-")
    # overall.S15.setextrainfo(["S15=25*S16"])
    # overall.S16 = overall.createstreams("S16")
    # overall.S16.setflow("-F")
    # overall.S16.setfractions([0.3,0.7])
    # # overall.S16.setequalto(splitter.S12)
    # print(overall)
    # overallSoln = overall.solvesystem()
    # print(overallSoln)
    
    # crystallizer.S7.setequalto(evaporator.S6, flowdirection="+")
    # # crystallizer.S7.setflow("F")
    
    # crystallizer.S8 = crystallizer.createstreams("S8")
    # crystallizer.S8.setequalto(overallSoln.streamregister["S15"], flowdirection="-")
    # crystallizer.S8.deleteextrainfo()
    # print(crystallizer)
    # crystallizerSoln = crystallizer.solvesystem()
    # print(crystallizerSoln)
    
    # splitter.S10.setequalto(crystallizerSoln.S9, flowdirection="+")
    
    # # splitter.S10.flowrate = -S10.flowrate
    # # splitter.S12.setequalto(overallSoln.streams["S16"])
    # splitter.S12.setequalto(overallSoln.S16, flowdirection="-")
    
    # print(splitter)
    # splitterSoln = splitter.solvesystem()
    # print(splitterSoln)
    
    # feedMixer.S2.setequalto(splitterSoln.S11, flowdirection="+")
    # # feedMixer.S2.flowrate = -S2.flowrate
    
    
    # print(feedMixer)
    # # print(feedMixer)
    # feedMixerSoln = feedMixer.solvesystem()
    # print(feedMixerSoln)
    # evaporator.S4.setequalto(feedMixerSoln.S3, flowdirection="+")
    # # evaporator.S4.flowrate = -S4.flowrate
    # print(evaporator)
    # evaporatorSoln = evaporator.solvesystem()
    # print(evaporatorSoln)
    
    # ##########################################
    # ##########################################
    
    # # usage for creating multiple streams using lists
    
    streams = ["Water", "NaOH", "Product"]
    flows = ["F", 1000, "-F"]
    fractions = [["x",0],
                  [0,1],
                  ["x", 0.1]]
    extrainfo = [[],[],["Water=0.9*Product"]]
    Mix = PhysicalProcess("A Mixer")
    # streamnames=[[]],flowrrates=[[]], fractions=[[]],extrainfos=[[]]
    
    Mix.createstreams(streamnames=streams, flowrrates=flows, fractions=fractions, extrainfos=extrainfo)
    # or use as follows. In the following case, each individual stream is captured in a variable
    # Mix.Water, Mix.NaOH, Mix.Product = Mix.createstreams(streamnames=streams, flowrrates=flows, fractions=fractions, extrainfos=extrainfo)

    
    print(Mix)
    print(Mix.solvesystem())
    Mix.degreesoffreedom()
    
    # usage for single stream
    streams = ["Water"]
    flows = ["F"]
    fractions = [["x",0]]
    extrainfo = [[]]
    Mix2 = PhysicalProcess("Another Mixer")
    Mix2.Water, = Mix2.createstreams(streamnames=streams, flowrrates=flows, fractions=fractions, extrainfos=extrainfo)
    print(Mix2)
    Mix2.HCl = Mix2.createstreams(streamname="HCl")
    Mix2.HCl.setequalto(-Mix2.Water)
    print(Mix2)
    
    
    # ##########################################
    # ##########################################
    # # Page 49, Example 2.5 Author: Sigurd Skogestad
    # # Create an instance of the physical process
    # # Use the following format
    # # variable name = pp.PhysicalProcess("Name of process")
    # mixer1 = PhysicalProcess("Ethanol Mixer")
    
    # # Start assigning streams 
    # mixer1.S1 = mixer1.createstreams("S1")
    # mixer1.S2 = mixer1.createstreams("S2")
    # mixer1.S3 = mixer1.createstreams("S3")
    
    # # Set data for streams
    # # Flowrate format: '+' if input stream and '-' if output stream
    # # Use the alphabet: "F" for input stream unknown flowrate
    # # Use the alphabet with negative added: "-F" for output stream unknown flowrate
    # # Do not add space between '-' and 'F', so "- F" is incorrect
    # # To give stream fractions use a list
    # # Use "x" to indicate an unknown fraction
    # mixer1.S1.setflow(-10)
    # mixer1.S1.setfractions([0.05, 1-0.05])
    # mixer1.S2.setflow("F")
    # mixer1.S2.setfractions([0.01, 1-0.01])
    # mixer1.S3.setflow("F")
    # mixer1.S3.setfractions([0.41, 1-0.41])
    # mixer1.sol = mixer1.solvesystem()
    # print(mixer1.sol)
    
    
    # ##########################################
    # ##########################################
    # # Page 49, Example 2.6 Author: Sigurd Skogestad
    # mixer2 = PhysicalProcess("NaCl mixer")
    # mixer2.S1 = mixer2.createstreams("S1")
    # mixer2.S2 = mixer2.createstreams("S2")
    # mixer2.S3 = mixer2.createstreams("S3")
    # mixer2.S1.setflow(-1)
    # mixer2.S1.setfractions([0.02, "x"])
    # mixer2.S2.setflow("F")
    # mixer2.S2.setfractions([0, 1])
    # mixer2.S3.setflow("F")
    # mixer2.S3.setfractions([.05, 1-0.05])
    # mixer2.sol = mixer2.solvesystem()
    # print(mixer2.sol)
    
    
    # ##########################################
    # ##########################################
    # # Page 51, Example 2.9 Author: Sigurd Skogestad
    # distColmn = PhysicalProcess("Distillation Column")
    # distColmn.S1 = distColmn.createstreams("S1")
    # distColmn.S2 = distColmn.createstreams("S2")
    # distColmn.S3 = distColmn.createstreams("S3")
    # distColmn.S4 = distColmn.createstreams("S4")
    # distColmn.S5 = distColmn.createstreams("S5")
    # distColmn.S1.setflow(800)
    # distColmn.S1.setfractions([0.7, 0.2933, 0.0067])
    # distColmn.S2.setflow("F")
    # distColmn.S2.setfractions([0.2, 0.8, 0])
    # distColmn.S3.setflow(-700)
    # distColmn.S3.setfractions([0.995, 0.005, 0])
    # distColmn.S4.setflow("-F")
    # distColmn.S4.setfractions([0.002, "x", 0.045,])
    # distColmn.S5.setflow("-F")
    # distColmn.S5.setfractions([0.002, 0.998, 0])
    # print(distColmn)
    # distColmn.sol = distColmn.solvesystem()
    # print(distColmn.sol)
    
    # ##########################################
    # ##########################################
    # column=PhysicalProcess("Liquid-Liquid Extraction and Distillation Column")
    
    # S1 = column.createstreams("S1")
    # S1.setflow("F")
    # #Component-1: Acetic Acid, component-2: Water, component-3: hexanol
    # S1.setfractions([0.18, 0.82, 0])
    
    # S2 = column.createstreams("S2") 
    # #Assuming flow rate m2 is 100
    # S2.setflow(100) 
    # S2.setfractions([0, 0, 1.00])
    # #S2.setextrainfo(["3:S2=(1/0.95)*S6"])
    
    # S4 = column.createstreams("S4")
    # S4.setflow("-F")
    # S4.setfractions([0.005, 0.995, 0])
    
    # S5 = column.createstreams("S5")
    # S5.setflow("-F")
    # S5.setfractions([0.96, 0, 0.04])
    
    # S6 = column.createstreams("S6")
    # S6.setflow("-F")
    # S6.setfractions([0.028, 0, 0.972])
    # S6.setextrainfo(["3:S6=0.95*S2"])
    
    # print(column)
    
    # ans = column.solvesystem()
    # print(ans)
    
