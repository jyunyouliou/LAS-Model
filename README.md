# LAS-Model

Focal seizure model
Author: Jyun-you Liou et al. 

To simulate focal seizures, simply run any of following files:

Exp1.m 
Exp2.m 
Exp4.m 
Exp5.m 
Exp6.m  
Exp7.m 

These scripts reproduce all the simulations we presented in the main text, corresponding with Figure 1, 2, 4, 5, 6, and 7, respectively.

For each experiment, the model's parameters are recorded in the corresponding template files, such as Exp1Template.m.

If you want to further modify the models, please read on.  Otherwise, the above is all you need to run the simulations. 

--------------- More information about the structure of this code --------------

%%%%%%%%%%%%%%% Class hierarchy %%%%%%%%%%%%%

A 'NeuralNetwork' can be classified into two major categories, 'RateNetwork' & 'SpikingNetwork'.  Both of which are sub-classes of NeuralNetwork

A 'NeuralNetwork' provides basic properties of a neural network, such as 'n' (number of neurons), 't' (current time), etc.  

Specific models are defined as sub-class of RateNetwork or SpikingNetwork.

Sub-classes of 'SpikingNetwork' are the following models:  
        % 'SpikingModel' - a biophysically motivated, spiking model.

Sub-classes of 'RateNetwork' are the following models: 
	% 'GeneralModel' - an abstract, rate based model. (Not used for this epilepsy study)
        % 'MeanFieldModel' - a mean-field approach of SpikingModel, using chloride concentration to determine inhibitory strength.
	% 'PhenomenologicalRateModel' - analogous to MeanFieldModel; however, local inhibition strength is described by an abstract variable 'z'. 

%%%%%%%% How NeuralNetwork computes %%%%%%%%%

Specifically how neurons inside a 'NeuralNetwork' compute, please look at each model's method 'IndividualModelUpdate'

All dynamical variables (those which change with time) are registered in the class definition file of each model.

A 'NeuralNetwork' can project/receive several 'Projection' to/from itself or other NeuralNetwork.  
Use NeuralNetwork's 'Link' method to create Projections (Or you can directly call 'Projection' from Class Projection. They are equivalent.)

The information is stored in .Proj.In for incoming projections and .Proj.Out for outward projections.
Proj.In & Proj.Out contains instances of 'Projection' Class  
	% 'Projection' contains the .Source (Pre-synaptic) and .Target (Post-synaptic) NeuralNetwork.
	% The properties of 'Projection' defines the type, topology, strength, and learning parameters ... of the synaptic connections  

Modify properties of 'Projection' instances to create various types of connections. (e.g. 'Topology','linear','Method','convolution')

%%%%%%%% How to start an experiment %%%%%%%%%

For one experiment, you need group all instances of NeuralNetwork into one MATLAB variable (vector) and use the method 'Update'.
  
'Update' will do the following things step by step:

1) Calculate the effective firing rate/spiking strength of NeuralNetwork and tell its .Proj.Out the result (stored in .Value of the 'Projection' instance)
2) The Projection instance will redistribute its effects (.Value) and tell the posy-synaptic NeuralNetwork what .Type of .Input they are receiving  

%%%%%%%%%%%%%% External Input %%%%%%%%%%%%%%%

You can create an external input by specify an instance of Class 'ExternalInput'.  

ExternalInput can be divided into two parts - .Deterministic and .Random.  See see help ExternalInput for more information how to set it up.

%%%%%%%%% Record Experiment Results %%%%%%%%%

A 'NeuralNetwork' can associate with a 'Recorder'.

	You can write dynamical variables & spikes into a 'Recorder' by NeuralNetwork's 'WriteToRecorder' method
	Depending on whether the model is a RateNetwork or SpikingNetwork, the recorder will automatically choose how to record it
	When the capacity of 'Recorder' is full, the program will detect it and automatically compress a 'Recorder' to another 'Recorder'
        Recorder.Ancestor indicates which NeuralNetwork their data originally come from, but Recorder.Parent tells which object 
        , either a NeuralNetwork or a Recorder, its data were directly derived from.  

%%%%%%% Plotting & Realtime feedback %%%%%%%%

You can directly plot an instance of 'NeuralNetwork' or 'Recorder' by their method 'plot'.

You can use 'AttachHotKey' (a method of NeuralNetwork) to your realtime simulation!  

	1) You need to plot the NeuralNetwork for hot keys to work. 
        2) After attach hot keys with the plot, you can give artificial stimulation, lock parameters, load previous data, restart, pause ... etc.

Here are a brief list of hot keys (it may change, please see help AttachHotKey to see detailed information 

Use NeuralNetwork's method 'AddHotKey' to allow real-time feedback from users to interact with the networks.  Here are the hot key lists:

'space' - pause/resume the simulation
'q' - quit the simulation
's' - give artificial external input as stimulation (left key is to give the waveform and right key is to give inversed polarity)
'c' - clamp dynamic variables (then left key is tuning up and right key is tuning down)
'p' - set parameter change rules
'uparrow' & 'downarrow' - increase or decrease parameters 

%%%%%%%%%%%%% Advanced Functions %%%%%%%%%%%%

For users who want to create their own models 

1) First, you need to decide whether it is a SpikingNetwork or a RateNetwork, and then inherit from it
   For SpikingNetwork, use the property .S (a structure by itself) to save spike & spike-derived variables
   For RateWork, use the property .R (also a structure by itself) to save rate & rate-derived variables

2) The method, 'IndividualModelUpdate' needs to be defined and requires the following operations
   a) Collect ExternalInput from .Ext & SynapticInput from Proj.In.Value
   b) Update internal dynamical variables
