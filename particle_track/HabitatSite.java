/*
 * Class to represent a site at which particles can be released/settle in particle tracking
 */
package particle_track;

//import java.util.Arrays;
import java.util.*;
import java.util.stream.IntStream;

/**
 *
 * @author SA01TA
 */
public class HabitatSite {
    private String ID; // A unique identifier for a site (string e.g. SEPA site IDs for fish farms
    private float[] xy; // Coordinate for site
    private float depth;
    private float scale; // A scale factor to be applied if required (e.g. site biomass for fish farms)
    private int containingMesh; // the mesh containing the habitat (possibly multiple in case of polygon?)
    private String containingMeshType;
    private int nearestFVCOMCentroid; // the mesh u,v point that is closest to this location
    private int containingFVCOMElem; // the element containing the habitat (possibly multiple in case of polygon?)
    private int[] nearestROMSGridPoint;
    private int[] containingROMSElem;
        
    public HabitatSite(String ID, float x, float y, float depth, float scale, List<Mesh> meshes)
    {
        this.ID = ID;
        this.xy = new float[]{x,y};
        this.depth = depth;
        this.scale = scale;
        
        double[] xy2 = new double[]{this.xy[0],this.xy[1]};
        
        // Find out which mesh the particle is in, based on bounding box.
        // We make the assumption that meshes appear in the list in their
        // order of precedence, and allocate particles to the first mesh that
        // contains their location
        this.containingMesh = -1;
        for (int m = 0; m < meshes.size(); m ++)
        {
            if (meshes.get(m).isInMesh(xy2))
            {
                this.containingMesh = m;
                break;
            }
        }
        if (this.containingMesh == -1)
        {
            System.err.println("Habitat site "+ID+" not within any provided mesh, check coordinates: "+x+" "+y);
        }
        
        this.nearestFVCOMCentroid = -1;
        this.containingFVCOMElem = -1;
        if (meshes.get(this.containingMesh).getType().equalsIgnoreCase("FVCOM"))
        {
            this.containingMeshType = "FVCOM";
            System.out.println("habitat site in FVCOM mesh");
            this.nearestFVCOMCentroid = Particle.nearestCentroid(xy[0], xy[1], meshes.get(this.containingMesh).getUvnode());
            this.containingFVCOMElem = Particle.whichElement(xy2, 
                        IntStream.rangeClosed(0, meshes.get(0).getUvnode()[0].length-1).toArray(), 
                        meshes.get(this.containingMesh).getNodexy(), meshes.get(this.containingMesh).getTrinodes());
        } 
        else if (meshes.get(this.containingMesh).getType().equalsIgnoreCase("ROMS"))
        {
            this.containingMeshType = "ROMS";
            System.out.println("habitat site in ROMS mesh");

            this.nearestROMSGridPoint = Particle.nearestROMSGridPoint(xy[0], xy[1], 
                    meshes.get(this.containingMesh).getLonU(), meshes.get(this.containingMesh).getLatU());
            
            System.out.println("found nearest point");
            // More to do here to turn the nearest grid point into the containing element
            this.containingROMSElem = Particle.whichROMSElement(xy[0], xy[1], 
                    meshes.get(this.containingMesh).getLonU(), meshes.get(this.containingMesh).getLatU(), 
                    this.nearestROMSGridPoint);
            System.out.println("found which element");
            
            
            
            
            
            
            System.out.println("Habitat site --- nearestROMSpoint=["+this.nearestROMSGridPoint[0]+","+this.nearestROMSGridPoint[1]+"]("
                    +meshes.get(this.containingMesh).getLonU()[this.nearestROMSGridPoint[0]][this.nearestROMSGridPoint[1]]+","
                    +meshes.get(this.containingMesh).getLatU()[this.nearestROMSGridPoint[0]][this.nearestROMSGridPoint[1]]+")"
                    +" --- containingROMSElem=["+this.containingROMSElem[0]+","+this.containingROMSElem[1]+"]");
        }
        
    }
    
    @Override
    public String toString()
    {
        String details = this.ID+" "+Arrays.toString(this.xy)+" "+this.containingMesh;
        if (this.containingMeshType.equalsIgnoreCase("FVCOM"))
        {
            details = this.ID+" "+Arrays.toString(this.xy)+" "+this.containingMesh+" "+this.nearestFVCOMCentroid+" "+this.containingFVCOMElem;
        }
        else if (this.containingMeshType.equalsIgnoreCase("ROMS"))
        {
            details = this.ID+" "+Arrays.toString(this.xy)+" "+this.containingMesh+" ("
                    +this.nearestROMSGridPoint[0]+","+this.nearestROMSGridPoint[1]+") ("
                    +this.containingROMSElem[0]+","+this.nearestROMSGridPoint[1]+")";
        }
        return details;
    }
    
    public String getID()
    {
        return this.ID;
    }
    
    public float[] getLocation()
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
    public int getContainingMesh()
    {
        return this.containingMesh;
    }
    public int getContainingFVCOMElem()
    {
        return this.containingFVCOMElem;
    }
    public int[] getContainingROMSElem()
    {
        return this.containingROMSElem;
    }

}
