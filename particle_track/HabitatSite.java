/*
 * Class to represent a site at which particles can be released/settle in particle tracking
 */
package particle_track;

import java.util.Arrays;

/**
 *
 * @author SA01TA
 */
public class HabitatSite {
    private String ID; // A unique identifier for a site (string e.g. SEPA site IDs for fish farms
    private double[] xy; // Coordinate for site
    private double depth;
    private double scale; // A scale factor to be applied if required (e.g. site biomass for fish farms)
    
    public HabitatSite(String ID, double x, double y, double depth, double scale)
    {
        this.ID = ID;
        this.xy = new double[]{x,y};
        this.depth = depth;
        this.scale = scale;
    }
    
    @Override
    public String toString()
    {
        return this.ID+" "+Arrays.toString(this.xy);
    }
    
    public String getID()
    {
        return this.ID;
    }
    
    public double[] getLocation()
    {
        return this.xy;
    }
    
    public double getDepth()
    {
        return this.depth;
    }
    
    public double getScale()
    {
        return this.scale;
    }

}
