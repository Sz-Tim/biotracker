/*
 * Class to represent a site at which particles can be released/settle in particle tracking
 */
package particle_track;

//import java.util.Arrays;
import java.util.*;
import java.util.stream.IntStream;
import static particle_track.Particle.distanceEuclid2;

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
    private boolean insideMesh;
    private String containingMeshType;
    private int nearestFVCOMCentroid; // the mesh u,v point that is closest to this location
    private int containingFVCOMElem; // the element containing the habitat (possibly multiple in case of polygon?)
    private int[] nearestROMSGridPointU;
    private int[] containingROMSElemU;
    private int[] nearestROMSGridPointV;
    private int[] containingROMSElemV;
        
    public HabitatSite(String ID, float x, float y, float depth, float scale, List<Mesh> meshes)
    {
        this.ID = ID;
        this.xy = new float[]{x,y};
        this.depth = depth;
        this.scale = scale;
        
        double[] xy2 = new double[]{this.xy[0],this.xy[1]};
        
        //System.out.println("### HABITAT SITE: "+ID+" ("+this.xy[0]+","+this.xy[1]+") ###");
        
        // Find out which mesh the particle is in, based on bounding box.
        // We make the assumption that meshes appear in the list in their
        // order of precedence, and allocate particles to the first mesh that
        // contains their location
        //this.containingMesh = -1;
        // Assume as a default that the habitat site is in the first mesh, if not change it to something else
        this.containingMesh = 0;
        this.insideMesh = false;
        
        //System.out.println("number of meshes "+meshes.size());
        for (int m = 0; m < meshes.size(); m ++)
        {
            //System.out.println("Mesh "+m+" "+meshes.get(m).isInMesh(xy2,true,new int[]{1}));
            if (meshes.get(m).isInMesh(xy2,true,null))
            {
                this.containingMesh = m;
                this.insideMesh = true;
                break;
            }
        }
        
        this.containingMeshType = "NONE";
        this.nearestFVCOMCentroid = -1;
        this.containingFVCOMElem = -1;
        
        //System.out.println("Containing mesh = "+this.containingMesh);
        if (this.insideMesh == false)
        {
            System.err.println("Habitat site "+ID+" not within any provided mesh --- defaulting to first mesh --- check coordinates: "+x+" "+y);
        } 
        else
        {
            if (meshes.get(this.containingMesh).getType().equalsIgnoreCase("FVCOM"))
            {
                this.containingMeshType = "FVCOM";
                //System.out.println("habitat site in FVCOM mesh");
                this.nearestFVCOMCentroid = Particle.nearestCentroid(xy[0], xy[1], meshes.get(this.containingMesh).getUvnode());
                if (this.insideMesh == false)
                {
                    double d1 = distanceEuclid2(xy[0], xy[1], meshes.get(this.containingMesh).getUvnode()[0][this.nearestFVCOMCentroid], meshes.get(this.containingMesh).getUvnode()[1][this.nearestFVCOMCentroid], "WGS84");
                    System.out.println("Habitat site ("+xy[0]+","+xy[1]+") outside mesh. Mearest centroid: "+this.nearestFVCOMCentroid
                                        +" ("+meshes.get(this.containingMesh).getUvnode()[0][this.nearestFVCOMCentroid]+","+meshes.get(this.containingMesh).getUvnode()[1][this.nearestFVCOMCentroid]+")"
                                        +" distance = "+d1);
                    if (d1 < 5000)
                    {
                        System.out.println("Within 5000m from mesh edge; moving to nearest element centroid");
                        this.xy = new float[]{meshes.get(this.containingMesh).getUvnode()[0][this.nearestFVCOMCentroid],meshes.get(this.containingMesh).getUvnode()[1][this.nearestFVCOMCentroid]};
                        xy2 = new double[]{this.xy[0],this.xy[1]};
                    }
                    else
                    {
                        System.err.println("More than 5000m from mesh edge; cannot create habitat site");
                    }                
                }
                this.containingFVCOMElem = Particle.whichElement(xy2, 
                            IntStream.rangeClosed(0, meshes.get(0).getUvnode()[0].length-1).toArray(), 
                            meshes.get(this.containingMesh).getNodexy(), meshes.get(this.containingMesh).getTrinodes());
            } 
            else if (meshes.get(this.containingMesh).getType().equalsIgnoreCase("ROMS"))
            {
                this.containingMeshType = "ROMS";
                System.out.println("habitat site in ROMS mesh");

                this.nearestROMSGridPointU = Particle.nearestROMSGridPoint(xy[0], xy[1], 
                        meshes.get(this.containingMesh).getLonU(), meshes.get(this.containingMesh).getLatU(), null);
                this.nearestROMSGridPointV = Particle.nearestROMSGridPoint(xy[0], xy[1], 
                        meshes.get(this.containingMesh).getLonV(), meshes.get(this.containingMesh).getLatV(), null);

                System.out.println("found nearest point");
                // More to do here to turn the nearest grid point into the containing element
                this.containingROMSElemU = Particle.whichROMSElement(xy[0], xy[1], 
                        meshes.get(this.containingMesh).getLonU(), meshes.get(this.containingMesh).getLatU(), 
                        this.nearestROMSGridPointU);
                this.containingROMSElemV = Particle.whichROMSElement(xy[0], xy[1], 
                        meshes.get(this.containingMesh).getLonV(), meshes.get(this.containingMesh).getLatV(), 
                        this.nearestROMSGridPointV);
                System.out.println("found which element");

                System.out.println("Habitat site --- nearestROMSpointU=["+this.nearestROMSGridPointU[0]+","+this.nearestROMSGridPointU[1]+"]("
                        +meshes.get(this.containingMesh).getLonU()[this.nearestROMSGridPointU[0]][this.nearestROMSGridPointU[1]]+","
                        +meshes.get(this.containingMesh).getLatU()[this.nearestROMSGridPointU[0]][this.nearestROMSGridPointU[1]]+")"
                        +" --- containingROMSElemU=["+this.containingROMSElemU[0]+","+this.containingROMSElemU[1]+"]");
                System.out.println("Habitat site --- nearestROMSpointV=["+this.nearestROMSGridPointV[0]+","+this.nearestROMSGridPointV[1]+"]("
                        +meshes.get(this.containingMesh).getLonV()[this.nearestROMSGridPointV[0]][this.nearestROMSGridPointV[1]]+","
                        +meshes.get(this.containingMesh).getLatV()[this.nearestROMSGridPointV[0]][this.nearestROMSGridPointV[1]]+")"
                        +" --- containingROMSElemV=["+this.containingROMSElemV[0]+","+this.containingROMSElemV[1]+"]");
            }
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
            details = this.ID+" "+Arrays.toString(this.xy)+" "+this.containingMesh+" U_grid: ("
                    +this.nearestROMSGridPointU[0]+","+this.nearestROMSGridPointU[1]+") ("
                    +this.containingROMSElemU[0]+","+this.nearestROMSGridPointU[1]+")";
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
    public String getContainingMeshType()
    {
        return this.containingMeshType;
    }
    public int getContainingFVCOMElem()
    {
        return this.containingFVCOMElem;
    }
    public int[] getContainingROMSElemU()
    {
        return this.containingROMSElemU;
    }
    public int[] getContainingROMSElemV()
    {
        return this.containingROMSElemV;
    }
    public int[] getNearestROMSPointU()
    {
        return this.nearestROMSGridPointU;
    }
    public int[] getNearestROMSPointV()
    {
        return this.nearestROMSGridPointV;
    }

}
