/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package tudelft.cgv.volume;

import java.io.File;
import java.io.IOException;

/**
 *
 * @author michel modified by Anna
 */

//////////////////////////////////////////////////////////////////////
///////////////// CONTAINS FUNCTIONS TO BE IMPLEMENTED ///////////////
//////////////////////////////////////////////////////////////////////

public class Volume {
    

    //Do NOT modify these attributes
    private int dimX, dimY, dimZ;
    private short[] data;
    private int[] histogram;

    // Do NOT modify this function
    // This function returns the nearest neighbour given a position in the volume given by coord.
    public short getVoxelNN(double[] coord) {
        if (coord[0] < 0 || coord[0] > (dimX-1) || coord[1] < 0 || coord[1] > (dimY-1)
                || coord[2] < 0 || coord[2] > (dimZ-1)) {
            return 0;
        }
        /* notice that in this framework we assume that the distance between neighbouring voxels is 1 in all directions*/
        int x = (int) Math.round(coord[0]); 
        int y = (int) Math.round(coord[1]);
        int z = (int) Math.round(coord[2]);
    
        return getVoxel(x, y, z);
    }
        

    //Do NOT modify this function
    //This function linearly interpolates the value g0 and g1 given the factor (t) 
    //the result is returned. It is used for the tri-linearly interpolation the values 
    private float interpolate(float g0, float g1, float factor) {
        float result = (1 - factor)*g0 + factor*g1;
        return result; 
    }
             
    //Do NOT modify this function
    // This function returns the trilinear interpolated value of the position given by  position coord.
    public float getVoxelLinearInterpolate(double[] coord) {
        if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
                || coord[2] < 0 || coord[2] > (dimZ-2)) {
            return 0;
        }
        /* notice that in this framework we assume that the distance between neighbouring voxels is 1 in all directions*/
        int x = (int) Math.floor(coord[0]); 
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);
        
        float fac_x = (float) coord[0] - x;
        float fac_y = (float) coord[1] - y;
        float fac_z = (float) coord[2] - z;

        float t0 = interpolate(getVoxel(x, y, z), getVoxel(x+1, y, z), fac_x);
        float t1 = interpolate(getVoxel(x, y+1, z), getVoxel(x+1, y+1, z), fac_x);
        float t2 = interpolate(getVoxel(x, y, z+1), getVoxel(x+1, y, z+1), fac_x);
        float t3 = interpolate(getVoxel(x, y+1, z+1), getVoxel(x+1, y+1, z+1), fac_x);
        float t4 = interpolate(t0, t1, fac_y);
        float t5 = interpolate(t2, t3, fac_y);
        float t6 = interpolate(t4, t5, fac_z);
        
        return t6; 
    }


    float a = -0.75f; // global variable that defines the value of a used in cubic interpolation.
    // you need to chose the right value
        
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
        
    // Function that computes the weights for one of the 4 samples involved in the 1D interpolation 
    public float weight (float x, Boolean one_two_sample)
    {
         float result=1.0f;
         
         // to be implemented
       
         return (float)result; 
   }
    
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
    // Function that computes the 1D cubic interpolation. g0,g1,g2,g3 contain the values of the voxels that we want to interpolate
    // factor contains the distance from the value g1 to the position we want to interpolate to.
    // We assume the out of bounce checks have been done earlier
    
    public float cubicinterpolate(float g0, float g1, float g2, float g3, float factor) {

        // Lets define our values of our kernel at points g0 to g1
        // kernel is centered in interpolation point  (g1 + factor)

        // distances from kernel origin to all samples
//        float k_dist_g0 = factor + 1;
//        float k_dist_g2 = 1 - factor;
//        float k_dist_g3 = 2 - factor;

        //  reduced forms of kernel values
        float c_0 = (float) (a* Math.pow(factor, 3)  -  2*a* Math.pow(factor, 2)  +  a*factor);
        float c_1 = (float) ((a+2) * Math.pow(factor, 3)  -  (a+3)* Math.pow(factor, 2)  +  1);
        float c_2 = (float) ((-a-2) * Math.pow(factor, 3)  +  (2*a+3)* Math.pow(factor, 2)  -  a*factor);
        float c_3 = (float) (-a* Math.pow(factor, 3)  +  a* Math.pow(factor, 2) );

        // value at interpolation point
        float result = c_0 * g0 + c_1 * g1 + c_2 * g2 + c_3 * g3;
                            
        return result; 
    }
        
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
    // 2D cubic interpolation implemented here. We do it for plane XY. Coord contains the position.
    // We assume the out of bounce checks have been done earlier
    public float bicubicinterpolateXY(double[] coord,int z) {


        /* notice that in this framework we assume that the distance between neighbouring voxels is 1 in all directions*/

        // get closest coordinates to our point
        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);

        // get our fractional distances from position in x and y directions
        float fac_x = (float) coord[0] - x;
        float fac_y = (float) coord[1] - y;

        // along x direction, construct 4 points interpolated x
        float bot_line_x = cubicinterpolate(getVoxel(x-1, y, z),
                                            getVoxel(x, y, z),
                                            getVoxel(x+1, y, z),
                                            getVoxel(x+2, y, z), fac_x);


        float top_line_x =  cubicinterpolate(getVoxel(x-1, y+1, z),
                                             getVoxel(x, y+1, z),
                                             getVoxel(x+1, y+1, z),
                                             getVoxel(x+2, y+1, z), fac_x);


        float before_bot_line_x = cubicinterpolate(getVoxel(x-1, y-1, z),
                                                   getVoxel(x, y-1, z),
                                                   getVoxel(x+1, y-1, z),
                                                   getVoxel(x+2, y-1, z), fac_x);


        float after_top_line_x =  cubicinterpolate(getVoxel(x-1, y+2, z),
                                                   getVoxel(x, y+2, z),
                                                   getVoxel(x+1, y+2, z),
                                                   getVoxel(x+2, y+2, z), fac_x);

        // along y direction, get final value
        float result =  cubicinterpolate(before_bot_line_x,
                                         bot_line_x,
                                         top_line_x,
                                         after_top_line_x, fac_y);

        return result; 

    }
            
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
    // 3D cubic interpolation implemented here given a position in the volume given by coord.
    
    public float getVoxelTriCubicInterpolate(double[] coord) {
        if (coord[0] < 1 || coord[0] > (dimX-3) || coord[1] < 1 || coord[1] > (dimY-3)
                || coord[2] < 1 || coord[2] > (dimZ-3)) {
            return 0;
        }


        /* notice that in this framework we assume that the distance between neighbouring voxels is 1 in all directions*/
        int z = (int) Math.floor(coord[2]);
        float fac_z = (float) coord[2] - z;


        // Get values for the 4 XY planes , along the z direction

        float before_bot_plane_z = bicubicinterpolateXY(coord, z-1);
        float bot_plane_z = bicubicinterpolateXY(coord, z);
        float top_plane_z = bicubicinterpolateXY(coord, z+1);
        float after_top_plane_z = bicubicinterpolateXY(coord,z+2 );
        // to be implemented              
        float result = cubicinterpolate(before_bot_plane_z,bot_plane_z,top_plane_z,after_top_plane_z, fac_z) ;
        result = (result > 0)?result: 0;
        return result; 
        

    }


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

	
    //Do NOT modify this function
    public Volume(int xd, int yd, int zd) {
        data = new short[xd*yd*zd];
        dimX = xd;
        dimY = yd;
        dimZ = zd;
    }
    //Do NOT modify this function
    public Volume(File file) {
        
        try {
            VolumeIO reader = new VolumeIO(file);
            dimX = reader.getXDim();
            dimY = reader.getYDim();
            dimZ = reader.getZDim();
            data = reader.getData().clone();
            computeHistogram();
        } catch (IOException ex) {
            System.out.println("IO exception");
        }
        
    }
    
    //Do NOT modify this function
    public short getVoxel(int x, int y, int z) {
    	int i = x + dimX*(y + dimY * z);
        return data[i];
    }
    
    //Do NOT modify this function
    public void setVoxel(int x, int y, int z, short value) {
    	int i = x + dimX*(y + dimY * z);
        data[i] = value;
    }
    
	//Do NOT modify this function
    public void setVoxel(int i, short value) {
        data[i] = value;
    }
    
    //Do NOT modify this function
    public short getVoxel(int i) {
        return data[i];
    }
    
	//Do NOT modify this function
    public int getDimX() {
        return dimX;
    }
    
    //Do NOT modify this function
    public int getDimY() {
        return dimY;
    }
    
    //Do NOT modify this function
    public int getDimZ() {
        return dimZ;
    }

    //Do NOT modify this function
    public short getMinimum() {
        short minimum = data[0];
        for (int i=0; i<data.length; i++) {
            minimum = data[i] < minimum ? data[i] : minimum;
        }
        return minimum;
    }
    
    //Do NOT modify this function
    public short getMaximum() {
        short maximum = data[0];
        for (int i=0; i<data.length; i++) {
            maximum = data[i] > maximum ? data[i] : maximum;
        }
        return maximum;
    }
 
    //Do NOT modify this function
    public int[] getHistogram() {
        return histogram;
    }
    
    //Do NOT modify this function
    private void computeHistogram() {
        histogram = new int[getMaximum() + 1];
        for (int i=0; i<data.length; i++) {
            histogram[data[i]]++;
        }
    }
}
