#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy

import matplotlib.pylab as plt
plt.switch_backend("Qt5Agg")
from SurfaceConic import SurfaceConic



class Beam2(Shadow.Beam):

    def __init__(self):
        self.beam = Shadow.Beam()


    def translation(self,qdist1):
        """
        :param qdist1: translation vector
        :return:
        """

        if numpy.array(qdist1).size != 3:
            raise Exception("Input must be a vector [x,y,z]")

        self.rays[0,:] += qdist1[0]
        self.rays[1,:] += qdist1[1]
        self.rays[2,:] += qdist1[2]


    def rotate(self,theta1,axis=1,rad=1):
        """
        :param theta1: the rotation angle in degrees (default=0)
        :param axis: The axis number (Shadow's column) for the rotation
                    (i.e, 1:x (default), 2:y, 3:z)
        :param file:
        :param rad: set this flag when theta1 is in radiants
        :return:
        """

        if not rad:
            theta1 = theta1 * numpy.pi / 180

        a1 = self.rays.copy()

        if axis == 1:
            torot = [2,3]
        elif axis == 2:
            torot = [1,3]
        elif axis == 3:
            torot = [1,2]


        costh = numpy.cos(theta1)
        sinth = numpy.sin(theta1)

        tstart = numpy.array([1,4,7,16])

        for i in range(len(tstart)):

            newaxis = axis + tstart[i] - 1
            newaxisi = newaxis - 1
            newtorot = torot + tstart[i] - 1
            newtoroti = newtorot -1

            self.rays[newtoroti[0],:] =  a1[newtoroti[0],:] * costh + a1[newtoroti[1],:] * sinth
            self.rays[newtoroti[1],:] = -a1[newtoroti[0],:] * sinth + a1[newtoroti[1],:] * costh
            self.rays[newaxisi]       =  a1[newaxisi,:]









    def traceGrating(self,oe,number,newcode=True): 
        if not newcode:
            return self.traceOE(oe,number)

        # ;
        # ; INPUTS
        # ;
        #
        fmirr         = oe.FMIRR    # 1
        p             = oe.T_SOURCE # 1000.0       # source-mirror
        q             = oe.T_IMAGE  # 300.0        # mirror-image
        alpha         = oe.ALPHA    # 0.0      # mirror orientation angle
        theta_grazing = (90.0-oe.T_INCIDENCE) * numpy.pi / 180  # 5e-3     # grazing angle, rad
        fcyl          = oe1.FCYL
        f_convex      = oe1.F_CONVEX
   
        print("fmirr = %s, p=%f, q=%f, alpha=%f, theta_grazing=%f rad, fcyl=%d"%\
              (fmirr,p,q,alpha,theta_grazing,fcyl))

        if fmirr == 2:
            ccc = SurfaceConic.initialize_as_ellipsoid_from_focal_distances(p,q,theta_grazing,cylindrical=fcyl,switch_convexity=f_convex)
            print(ccc)
        else:
            raise Exception("Not implemented")



        #
        # put beam in mirror reference system
        #
        # TODO: calculate rotation matrices? Invert them for putting back to the lab system?
    
        self.rotate(alpha,axis=2)
        self.rotate(theta_grazing,axis=1)
        self.translation([0.0,-p*numpy.cos(theta_grazing),p*numpy.sin(theta_grazing)])
    
    
        #
        # reflect beam in the mirror surface and dump mirr.01
        #
        self = ccc.apply_specular_reflection_on_beam(self)
        #self.dump_shadow3_file('minimirr.02')
    
        #
        # put beam in lab frame and compute image
        #
        #self.rotate(theta_grazing,axis=1)
        # TODO what about alpha?
        #self.retrace(q,resetY=True)
        #self.dump_shadow3_file('ministar.02')















        return self.traceOE(oe,number)
    

	
if __name__ == "__main__":
    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0
    
    #
    # initialize shadow3 source (oe0) and beam
    #
    ##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    beam = Beam2()
    ##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    oe0 = Shadow.Source()
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()
    
    #
    # Define variables. See meaning of variables in: 
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml 
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #
    
    oe0.FDISTR = 3
    oe0.F_COLOR = 2
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.0
    oe0.HDIV2 = 0.0
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.NPOINT = 100000
    oe0.N_COLOR = 1
    oe0.PH1 = 1000.0
    oe0.PH2 = 1000.2
    oe0.SIGDIZ = 0.001
    oe0.SIGMAX = 0.0015
    oe0.SIGMAZ = 0.0015
    oe0.VDIV1 = 0.0
    oe0.VDIV2 = 0.0
    
    oe1.DUMMY = 0.1
    oe1.FWRITE = 3
    oe1.F_REFRAC = 2
    oe1.F_SCREEN = 1
    oe1.I_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe1.N_SCREEN = 1
    oe1.RX_SLIT = numpy.array([20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe1.RZ_SLIT = numpy.array([4.77, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe1.T_IMAGE = 0.0
    oe1.T_INCIDENCE = 0.0
    oe1.T_REFLECTION = 180.0
    oe1.T_SOURCE = 4300.0
    
    oe2.DUMMY = 0.1
    oe2.FCYL = 1
    oe2.FMIRR = 2
    oe2.FWRITE = 1
    oe2.F_DEFAULT = 0
    oe2.F_GRATING = 1
    oe2.F_RULING = 5
    oe2.RULING = 3500.0
    oe2.RUL_A1 = -3.1
    oe2.RUL_A2 = 0.00215
    oe2.SIMAG = 2800.0
    oe2.SSOUR = 4300.0
    oe2.THETA = 89.484
    oe2.T_IMAGE = 2136.0
    oe2.T_INCIDENCE = 89.484
    oe2.T_REFLECTION = 84.635452
    oe2.T_SOURCE = 0.0
    
    
    
    #Run SHADOW to create the source
    
    if iwrite:
        oe0.write("start.00")
    
    beam.genSource(oe0)
    
    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")
    
    
    #
    #run optical element 1
    #
    print("    Running optical element: %d"%(1))
    if iwrite:
        oe1.write("start.01")
    beam.traceOE(oe1,1)
    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")
    
    
    #
    #run optical element 2
    #
    print("    Running optical element: %d"%(2))
    if iwrite:
        oe2.write("start.02")

    ##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #beam.traceOE(oe2,2)
    beam.traceGrating(oe2,2)
    ##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if iwrite:
        oe2.write("end.02")
        beam.write("star.02")
    
    
    Shadow.ShadowTools.plotxy(beam,1,3,nbins=101,nolost=1,title="Real space",xrange=[-15,15],yrange=[-0.010,0.010])
        
