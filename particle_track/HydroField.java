/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package particle_track;

import java.io.IOException;
import ucar.ma2.InvalidRangeException;

/**
 *
 * @author SA01TA
 */
public class HydroField {
    private float[][][] u;
    private float[][][] v;
    private float[][][] s;
    private float[][][] t;
    private float[][][] el;
    
    /**
     * Default constructor for a single day's data
     * 
     * @param filename 
     */
    public HydroField(String filename)
    {
        try
        {
            u = IOUtils.readNetcdfFloat3D(filename,"u");
            v = IOUtils.readNetcdfFloat3D(filename,"v");
            s = IOUtils.readNetcdfFloat3D(filename,"s");
            t = IOUtils.readNetcdfFloat3D(filename,"t");
            el = IOUtils.readNetcdfFloat3D(filename,"el");
            
        } catch (IOException ioe) {
	    ioe.printStackTrace();
        } catch (InvalidRangeException ire) {
            ire.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }   
    }
    
}
