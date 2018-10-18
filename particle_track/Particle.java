/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package particle_track;

import java.util.*;

/**
 *
 * @author tomdude
 */
public class Particle {
    
    final private int id;
    // horizontal position
    private double[] xy = new double[2];
    final private double[] startLoc = new double[2];
    private String startSiteID = "0";
    
    private ISO_datestr startDate;
    private double startTime = 0;
    
    private int elem;
    private double[][] nrList = new double[5][2];
    private double[][] cornerList = new double[3][2];
    // release time
    private double releaseTime = 0;
    // vertical position
    private double z = 0;
    private int depLayer = 0;
    // settlement details
    private double age = 0;
    private int status = 0;
    private double density = 1;
    private double mortRate = 0.01; // default hourly rate, based on results of Stein et al. 2005 (copepodid, nauplii rate marginally lower = 0.0078)
    private boolean arrived = false;
    private boolean viable = false;
    private boolean free = false;
    private boolean settledThisHour = false;
    private boolean boundaryExit = false;
    
    private String lastArrival = "0";
    
    // A list to store data on the arrivals made by each particle.
    // If rp.endOnArrival=true, this will contain a maximum of one element.
    // --- Not presently used ---
    private List<Arrival> arrivals;
    
    // create a new particle at a defined location, at the water surface
    public Particle(double xstart, double ystart, String startSiteID, int id, double mortalityRate)
    {
        this.id = id;
        this.xy[0] = xstart;
        this.xy[1] = ystart;
        this.startSiteID = startSiteID;
        this.startLoc[0] = xstart;
        this.startLoc[1] = ystart;
        this.mortRate = mortalityRate;
        
        this.arrivals = new ArrayList<>();
    }
    public Particle(double xstart, double ystart, String startSiteID, int id, double mortalityRate, 
            ISO_datestr startDate, double startTime)
    {
        this.id = id;
        this.xy[0] = xstart;
        this.xy[1] = ystart;
        this.startSiteID = startSiteID;
        this.startDate = startDate;
        this.startTime = startTime;
        this.startLoc[0] = xstart;
        this.startLoc[1] = ystart;
        this.mortRate = mortalityRate;
        
        //this.arrivals = new ArrayList<Arrival>();
    }

//    public void setReleaseScenario(int releaseScenario, double[][] startlocs)
//    {
//        switch (releaseScenario) {
//            // integer to switch release scenario
//            // 0 all at time zero
//            // 1 tidal release (evenly over first 24 hours)
//            // 2 continuous release (1 per hour per site)
//            // 3 continuous release (5 per hour per site)
//            // 4 continuous release (10 per hour per site)
//            // 5 continuous release (20 per hour per site)
//            // 10 defined release times
//            case 0:
//                this.setReleaseTime(0);
//                break;
//            case 1:
//                this.setReleaseTime((this.id / startlocs.length) % 25);
//                break;
//            case 2:
//                this.setReleaseTime(Math.floor(this.id / startlocs.length));
//                break;
//            case 3:
//                this.setReleaseTime(Math.floor(this.id / (5 * startlocs.length)));
//                break;
//            case 4:
//                this.setReleaseTime(Math.floor(this.id / (10 * startlocs.length)));
//                break;
//            case 5:
//                this.setReleaseTime(Math.floor(this.id / (20 * startlocs.length)));
//                break;
//            case 10:
//                this.setReleaseTime(startlocs[this.id][3]);
//                break;
//        }
//    }
    
    @Override
    public String toString()
    {
        return this.getID()+" "+this.xy.toString();
    }
    
    // Not presently used
    public void reportArrival(int sourceLocation, int arrivalLocation, double arrivalTime, double arrivalDensity)
    {
        arrivals.add(new Arrival(sourceLocation,arrivalLocation,arrivalTime,arrivalDensity));
        //System.out.printf("Arrival (particle %d): %d %d %f %f\n", this.getID(),sourceLocation,arrivalLocation,arrivalTime,arrivalDensity);
    }
    public List<Arrival> getArrivals()
    {
        return this.arrivals;
    }
    
    public void setReleaseTime(double releaseTime)
    {
        this.releaseTime = releaseTime;
    }
    public double getReleaseTime()
    {
        return this.releaseTime;
    }
    public int getID()
    {
        return this.id;
    }
    public double[] getLocation()
    {
        return this.xy;
    }
    public double[] getStartLocation()
    {
        return this.startLoc;
    }
    public String getStartID()
    {
        return this.startSiteID;
    }
    public ISO_datestr getStartDate()
    {
        return this.startDate;
    }
    public double getStartTime()
    {
        return this.startTime;
    }
    public void setLocation(double x, double y)
    {
        this.xy[0] = x;
        this.xy[1] = y;
    }
    public void setElem(int elem)
    {
        this.elem = elem;
    }
    public int getElem()
    {
        return this.elem;
    }
    public void setNrList(double[][] uvnode)
    {
        this.nrList = nearestCentroidList(this.xy[0], this.xy[1], uvnode);
    }
    public double[][] getNrList()
    {
        return this.nrList;
    }
    public void setCornerList(int elemPart, int[][] trinodes, double[][] nodexy)
    {
        this.cornerList[0][0]=nodexy[trinodes[elemPart][0]][0];
        this.cornerList[0][1]=nodexy[trinodes[elemPart][0]][1];
        this.cornerList[1][0]=nodexy[trinodes[elemPart][1]][0];
        this.cornerList[1][1]=nodexy[trinodes[elemPart][1]][1];
        this.cornerList[2][0]=nodexy[trinodes[elemPart][2]][0];
        this.cornerList[2][1]=nodexy[trinodes[elemPart][2]][1];
    }
    public double[][] getCornerList()
    {
        return this.cornerList;
    }
    
    public double getZ()
    {
        return this.z;
    }
    public void setZ(double z)
    {
        this.z = z;
    }
    public int getDepthLayer()
    {
        return this.depLayer;
    }
    public void setDepthLayer(int dep)
    {
        this.depLayer = dep;
    }
    // see lower down for setMortRate (based on salinity)
    public double getMortRate()
    {
        return mortRate;
    }
    public void setDensity()
    {
        this.density = this.density*(1-this.mortRate);
        //System.out.println("density = "+this.density);
    }
    public double getDensity()
    {
        return this.density;
    }
    public void setStatus(int status)
    {
        // 0 = not free
        // 1 = free
        // 2 = viable (able to settle)
        // 3 = settled
        // 66 = boundary exit
        // 666 = dead (exceeded lifespan)
        this.status=status;
    }
    public int getStatus()
    {
        return this.status;
    }
    
    /** Set particle depth in the water column based on its defined behaviour and
     *  the time
     * 1 - passive, stay on surface
     * 2 - passive, stay on bottom (layer 10)
     * 3 - passive, stay in mid layer (layer 5)
     * 4 - vertical swimming: surface for hours 19-6, mid layer (5) hours 7-18
     * 5 - rapid drop (1->10) at hour 6, then gradually move back up
     * 6 - top during flood tides, mid during ebb (local)
     * 7 - mid during flood tides, bed during ebb (local)
     * 8 - top during flood tides, bed during ebb (local)
     * MORE...? "homing" ability
     */
    public void setDepthLayer(int behaviour, String tideState)
    {
        boolean reducedFile = true;
        switch (behaviour)
        {
            case 1: this.depLayer = 0; break;
            case 2: this.depLayer = 9; break;
            case 3: this.depLayer = 5; break;
            case 6: 
                if (tideState.equalsIgnoreCase("flood"))
                {
                    this.depLayer = 0;
                } else
                {
                    this.depLayer = 5;
                }
                break;
            case 7: 
                if (tideState.equalsIgnoreCase("flood"))
                {
                    this.depLayer = 5;
                } else
                {
                    this.depLayer = 9;
                }
                break;
            case 8: 
                if (tideState.equalsIgnoreCase("flood"))
                {
                    this.depLayer = 0;
                } else
                {
                    this.depLayer = 9;
                }
                break;
        }
        // enable running file with two depth layers ("top"=0 and "bottom"=1)
        if (reducedFile==true)
        {
            if (this.depLayer>0)
            {
                this.depLayer=1;
            }
        }
        //return 0;
    }
   
    /** put particle in the correct depth layer, based upon
     * its z position, and
     * 
     * @param localDepth  element bathymetry value
     * @param layers 
     */
    public void setLayerFromDepth(double localDepth, double[] layers)
    {
        int depNearest = 0;
        double dZmin = 1000;
        for (int i = 0; i < layers.length; i++)
        {
            if (Math.abs(this.z - localDepth*layers[i]) < dZmin)
            {
                depNearest = i;
            }
        }
        //System.out.printf("setting depth layer: %d (%f, particle depth = %f)\n",depNearest,localDepth*layers[depNearest],this.z);
        this.setDepthLayer(depNearest);
    }
    
    public double verticalDiffusion()
    {
        double vertDiff = 0;
        return vertDiff;
    }
    
    
    public void setViable(boolean viable)
    {
        this.viable = viable;
    }
    public void setArrived(boolean arrived)
    {
        this.arrived = arrived;
    }
    public void setLastArrival(String loc)
    {
        this.lastArrival = loc;
    }
    public String getLastArrival()
    {
        return this.lastArrival;
    }
    public void setSettledThisHour(boolean settled)
    {
        this.settledThisHour = settled;
    }
    public void setFree(boolean free)
    {
        this.free = free;
    } 
    public void setBoundaryExit(boolean exit)
    {
        this.boundaryExit = exit;
    }
    public boolean getViable()
    {
        return this.viable;
    }
    public boolean getArrived()
    {
        return this.arrived;
    }
    public boolean getSettledThisHour()
    {
        return this.settledThisHour;
    }
    public boolean getFree()
    {
        return this.free;
    }
    public boolean getBoundaryExit()
    {
        return this.boundaryExit;
    }    
    public double[] behaveVelocity(int behaviour)
    {
        double[] uv = new double[2];
        return uv;
    }
    
    public double[] smagorinskyDiffusionVelocity(int node, int[] neighbours, double u, double v, double[][] uvnode)
    {
        double[] uv = new double[2];
        return uv;
    }
    
    /**
     * Compute salinity from an average of the corner nodes of the containing element.
     * @param tt
     * @param salinity
     * @param trinodes
     * @return 
     */
    public double salinity(int tt, double[][] salinity, int[][] trinodes)
    {
        double s = 0;
        double dist1 = distanceEuclid(this.xy[0],this.xy[1],cornerList[0][0],cornerList[0][1]);
        double dist2 = distanceEuclid(this.xy[0],this.xy[1],cornerList[1][0],cornerList[1][1]);
        double dist3 = distanceEuclid(this.xy[0],this.xy[1],cornerList[2][0],cornerList[2][1]);
        double weight1 = 1.0/(dist1*dist1);
        double weight2 = 1.0/(dist2*dist2);
        double weight3 = 1.0/(dist3*dist3);
        double weightSum = weight1+weight2+weight3;
        
        s = (1.0/weightSum)*(weight1*salinity[tt][trinodes[elem][0]]+weight2*salinity[tt][trinodes[elem][1]]+weight3*salinity[tt][trinodes[elem][2]]);
        
        return s;
    }
    /**
     * Sets the mortality rate for the particle packet based upon local salinity
     * @param salinity
     * @return 
     */
    public void setMortRate(double salinity)
    {
        // estimated 2nd order polynomial fit to Bricknell et al.'s (2006) data
        this.mortRate = 0.0011*salinity*salinity - 0.07*salinity + 1.1439;
        //System.out.println("salinity = "+salinity+" mortrate = "+this.mortRate);
    }
     
    public void increaseAge(double increment)
    {
        this.age+=increment;
    }
    public double getAge()
    {
        return this.age;
    }
    
    
    
    // --------------------------------------------------------------------------------------------------------------
    // Everything to do with velocity calculation below here
    // --------------------------------------------------------------------------------------------------------------
    /**
     * Find the nearest mesh element centroid
     * @param x
     * @param y
     * @param uvnode
     * @return 
     */
    public static int nearestCentroid(double x, double y, double[][] uvnode)
    {
        int nearest = -1;
        double dist=10000000;
    
        for (int i = 0; i < uvnode.length; i++)
        {
            double distnew = Math.sqrt((x-uvnode[i][0])*(x-uvnode[i][0])+(y-uvnode[i][1])*(y-uvnode[i][1]));
            if (distnew<dist)
            {
                dist=distnew;
                nearest=i;
            }
        }
        //System.out.printf("In Particle.nearestCentroid "+nearest+"\n");
        return nearest;
    }
    /**
     * Make a list of nearest mesh element centroids
     * @param x
     * @param y
     * @param uvnode
     * @return 
     */
    public static double[][] nearestCentroidList(double x, double y, double[][] uvnode)
    {
        double[][] nearestList = new double[5][2];
        int nearest = -1;
        double dist=10000000;
    
        for (int i = 0; i < uvnode.length; i++)
        {
            double distnew = Math.sqrt((x-uvnode[i][0])*(x-uvnode[i][0])+(y-uvnode[i][1])*(y-uvnode[i][1]));
            if (distnew<dist)
            {
                dist=distnew;
                nearest=i;
                // Shift everything along one element
                nearestList[4][0]=nearestList[3][0];
                nearestList[4][1]=nearestList[3][1];
                nearestList[3][0]=nearestList[2][0];
                nearestList[3][1]=nearestList[2][1];
                nearestList[2][0]=nearestList[1][0];
                nearestList[2][1]=nearestList[1][1];             
                nearestList[1][0]=nearestList[0][0];
                nearestList[1][1]=nearestList[0][1];

                nearestList[0][0]=i;
                nearestList[0][1]=dist;              
            }
        }
        return nearestList;
    }
    
    /**
     * Calculate velocity at a particle's location, given known
     * @param tt
     * @param u
     * @param v
     * @return 
     */
    public double[] velPart(int tt, double u[][], double v[][])
    {
        double[] velocity = new double[2];
        velocity = velocityFromNearestList(this.nrList,tt,u,v,10,1);
        return velocity;
    }
    /**
     * Method to compute velocity at a certain time and location, directly from an element of a velocity field.
     * 
     * This does the same thing as "velPart", except it generates a list of cells
     * to spatially interpolate from.
     * 
     * @param tt
     * @param xy
     * @param elemPart0
     * @param u
     * @param v
     * @param neighbours
     * @param uvnode
     * @param nodexy
     * @param trinodes
     * @param allelems
     * @return 
     */
    public static double[] velocityInterpAtLoc(int tt, double xy[], int elemPart0, double u[][],double v[][],
            int[][] neighbours, double[][] uvnode, double[][] nodexy, int[][] trinodes, int[] allelems, int depLayer)
    {
        double[] velocity = new double[2];
        // elemPart0 is the starting search element. Set = 1 if outside range.
        if (elemPart0 < 1 || elemPart0 > allelems.length)
        {
            elemPart0 = 1;
            System.out.println("velocityInterpAtLoc picking default start element = 1");
        }
        double[][] cellList = neighbourCellsList(xy, elemPart0, neighbours, uvnode, nodexy, trinodes, allelems);
        velocity = velocityFromNearestList(cellList,tt,u,v,10,1);
        return velocity;
    }
    
    /**
     * Underlying static method used by both the above two methods
     * @param nrList
     * @param tt
     * @param u
     * @param v
     * @param numLayers
     * @param depLayer
     * @return 
     */
    public static double[] velocityFromNearestList(double[][] nrList, int tt, double u[][],double v[][], int numLayers, int depLayer)
    {
        double[] velocity = new double[2];
        double[] weights = new double[nrList.length];
        double usum=0,vsum=0,sum=0;
        for (int i = 0; i < nrList.length; i++)
        {
            if (nrList[i][1]!=0)
            {
                weights[i]=1.0/(nrList[i][1]*nrList[i][1]);
            }
            else
            {
                weights[i]=1;
            }
            
            boolean reducedFile = true;
            if (reducedFile==true)
            {
                numLayers = 1;
            }
            //System.out.printf("tt %d i %d",tt,i);
            //System.out.printf(" --- elem %d dist %.4f weight %.4f --- vel = %.4f %.4f\n",(int)this.nrList[i][0],this.nrList[i][1],weights[i],u[tt][(int)this.nrList[i][0]*numLayers+this.depLayer],v[tt][(int)this.nrList[i][0]*numLayers+this.depLayer]);
            usum=usum+weights[i]*u[tt][(int)nrList[i][0]*numLayers+depLayer];
            vsum=vsum+weights[i]*v[tt][(int)nrList[i][0]*numLayers+depLayer];
            sum=sum+weights[i];
        }
        usum=usum/sum;
        vsum=vsum/sum;
        velocity[0]=usum;
        velocity[1]=vsum;
        //System.out.printf("Interpolated Velocity = %.4f %.4f\n",velocity[0],velocity[1]);
        return velocity;
    }
    /**
     * Find which element a particle resides within (edge checking)
     * @param x
     * @param y
     * @param elems
     * @param nodexy
     * @param trinodes
     * @return 
     */
    public static int whichElement(double x, double y, int[] elems, double[][] nodexy, int[][] trinodes)
    {
        int which = -1;
        int res = 0;
        for (int i = 0; i < elems.length; i++)
        {
            double[] xt = new double[3];
            double[] yt = new double[3];
            
            for (int j = 0; j < 3; j++)
            {
                //int elem=elems[i];
                
                 xt[j]=nodexy[trinodes[elems[i]][j]][0];
                 yt[j]=nodexy[trinodes[elems[i]][j]][1];
            }
            // check whether (x,y) lies within this
            //fprintf('check %d\n', possibleElems(i));
            
            double f1 = (y-yt[0])*(xt[1]-xt[0]) - (x-xt[0])*(yt[1]-yt[0]);
            double f2 = (y-yt[2])*(xt[0]-xt[2]) - (x-xt[2])*(yt[0]-yt[2]);
            double f3 = (y-yt[1])*(xt[2]-xt[1]) - (x-xt[1])*(yt[2]-yt[1]);
            if(f1*f3 >= 0.0 && f3*f2 >= 0.0) 
            {
                res = 1;
            }
            if(res==1)
            {
                which=elems[i];
                
                break;
            }

        }
        //System.out.printf("whichElement: particle in %d\n", which);
        return which;
    }
    /**
     * Euclidean distance between two points
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @return 
     */
    public static double distanceEuclid(double x1, double y1, double x2, double y2)
    {
        double dist = Math.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
        return dist;
    }
    /**
     * Search through lists of elements progressively further away from the last known
     * location in order to find the new location.
     * 
     * Returned array contains the element in which the particle is determined to be located, 
     * plus the number of counts required at each scale (for diagnostics).
     * @param newlocx
     * @param newlocy
     * @param elemPart
     * @param nodexy
     * @param trinodes
     * @param neighbours
     * @param allelems
     * @return 
     */
    public static int[] findContainingElement(double newlocx, double newlocy, int elemPart,
            double[][] nodexy, int[][] trinodes, int[][] neighbours, int[] allelems)
    {
        int[] c = new int[6];
        int[] elems = new int[1];
        elems[0] = elemPart;
        //System.out.println("findContainingElement: elems[0]="+elems[0]);
        int whereami=Particle.whichElement(newlocx,newlocy,elems,nodexy,trinodes);
        //System.out.println("findContainingElement: whereami="+whereami);
        c[1]=1;
        if (whereami==-1)
        {
            int[] elems0 = neighbours[elemPart];
            c[2]=1;
            whereami=Particle.whichElement(newlocx,newlocy,elems0,nodexy,trinodes);
            // if fails, look in nearest 10 (id numerical)
            if (whereami==-1)
            {
                c[3]=1;
                int[] elems1 = new int[10];
                for (int j = 0; j < 10; j++)
                {
                    elems1[j] = Math.min(Math.max(elemPart-5+j,0),allelems.length-1);
                }
                whereami=Particle.whichElement(newlocx,newlocy,elems1,nodexy,trinodes);
                // if fails, look in nearest 500 (id numerical)
                if (whereami==-1)
                {
                    c[4]=1;
                    int[] elems2 = new int[500];
                    for (int j = 0; j < 500; j++)
                    {
                        elems2[j] = Math.min(Math.max(elemPart-250+j,0),allelems.length-1);
                    }
                    whereami=Particle.whichElement(newlocx,newlocy,elems2,nodexy,trinodes);
                    // if this fails, look in all elements
                    if (whereami==-1)
                    {
                        c[5]=1;
                        whereami=Particle.whichElement(newlocx,newlocy,allelems,nodexy,trinodes);
                    }
                }
            }
        }
        c[0]=whereami;
        //System.out.printf("%d %d %d %d %d %d\n",c[0],c[1],c[2],c[3],c[4],c[5]);
//        if (c[0] == 0 || c[0] == -1)
//        {
//            System.out.println("Element out of bounds, fixing location");
//            System.out.printf("whereami=0 --- %.6e %.6e %d\n",newlocx,newlocy,elemPart);
//        }
        return c;
    }
    
     
    /**
     * Set the particle's nearestList (NrList) to be identical to the list of elements which 
     * neighbour the element that the particle is actually in
     * @param neighbours
     * @param uvnode 
     */
    public void setNrListToNeighbourCells(int[][] neighbours, double[][] uvnode)
    {      
        // distance to elem
        //int elem = nearestCentroid(this.xy[0],this.xy[1],uvnode);
        int elem = this.elem;
        this.nrList[0][0] = elem;
        this.nrList[0][1] = distanceEuclid(this.xy[0],this.xy[1],uvnode[elem][0],uvnode[elem][1]);
        // distance to neighbouring elems
        this.nrList[1][0] = neighbours[elem][0];
        this.nrList[1][1] = distanceEuclid(this.xy[0],this.xy[1],uvnode[neighbours[elem][0]][0],uvnode[neighbours[elem][0]][1]);
        this.nrList[2][0] = neighbours[elem][1];
        this.nrList[2][1] = distanceEuclid(this.xy[0],this.xy[1],uvnode[neighbours[elem][1]][0],uvnode[neighbours[elem][1]][1]);
        this.nrList[3][0] = neighbours[elem][2];
        this.nrList[3][1] = distanceEuclid(this.xy[0],this.xy[1],uvnode[neighbours[elem][2]][0],uvnode[neighbours[elem][2]][1]);   
        this.nrList[4][0] = 0;
        this.nrList[4][1] = 1000000;     
    }
    
    public static double[][] neighbourCellsList(double xy[], int elemPart0, 
            int[][] neighbours, double[][] uvnode, double[][] nodexy, int[][] trinodes, int[] allelems)
    {      
        double[][] nrList = new double[5][2];
        // distance to elem
        //int elem = nearestCentroid(this.xy[0],this.xy[1],uvnode);
        int elem[] = findContainingElement(xy[0], xy[1], elemPart0,
            nodexy, trinodes, neighbours, allelems);
        // If particle is not within the mesh (value returned by findContainingElement = -1)
        // exit this method returning array of zeros.
        if (elem[0] == -1)
        {
            return nrList;
        }
        int thisElem = elem[0]; 
        //int elem = this.elem;
        nrList[0][0] = thisElem;
        nrList[0][1] = distanceEuclid(xy[0],xy[1],uvnode[thisElem][0],uvnode[thisElem][1]);
        // distance to neighbouring elems
        nrList[1][0] = neighbours[thisElem][0];
        nrList[1][1] = distanceEuclid(xy[0],xy[1],uvnode[neighbours[thisElem][0]][0],uvnode[neighbours[thisElem][0]][1]);   
        nrList[2][0] = neighbours[thisElem][1];
        nrList[2][1] = distanceEuclid(xy[0],xy[1],uvnode[neighbours[thisElem][1]][0],uvnode[neighbours[thisElem][1]][1]);    
        nrList[3][0] = neighbours[thisElem][2];
        nrList[3][1] = distanceEuclid(xy[0],xy[1],uvnode[neighbours[thisElem][2]][0],uvnode[neighbours[thisElem][2]][1]);      
        nrList[4][0] = 0;
        nrList[4][1] = 1000000; 
        
        return nrList;
    }
    
    /** Method for updating location using an RK4 integration step
     * 
     * Note that this must use a spatially interpolated velocity for each step 
     * of the calculation. Therefore read in full velocity fields here
     * 
     * 14/02/17 --- Removed  double[][] u1, double[][] v1, from argument list. 
     *              Removed from method calls stepAhead
     *              Replace velplus1[] calculation with direct calculation
     */
    public double[] rk4Step(double[][] u, double[][] v, // velocities
            int[][] neighbours, double[][] uvnode, double[][] nodexy, 
            int[][] trinodes, int[] allelems,      // other mesh info
            int tt, int st, double dt,                                  // locate particle in space and time
            int stepsPerStep, int numLayers)   // info on simulation length
    {
        int elemPart = this.getElem();
        int dep = this.getDepthLayer();
        //System.out.printf("RK4Step: Location = [%.6e,%.6e], Element = %d\n",this.getLocation()[0],this.getLocation()[1],elemPart);
        double[] advectStep = new double[2];
        this.setNrListToNeighbourCells(neighbours,uvnode);
        
        // 2. Compute k_1 (spatial interpolation at start of step)
        //System.out.println("Start step");
        double[] vel = this.velPart(tt,u,v);
        double[] velplus1 = velocityFromNearestList(this.getNrList(),tt+1,u,v,numLayers,dep);
        //double[] velplus1 = getNextVel(tt,recordsPerFile,fnum,lastday,this.getNrList(),u,v,u1,v1,numLayers,dep);
        double[] k1 = new double[2];
        k1[0] = dt*(vel[0] + ((double)st/(double)stepsPerStep)*(velplus1[0]-vel[0]));
        k1[1] = dt*(vel[1] + ((double)st/(double)stepsPerStep)*(velplus1[1]-vel[1]));
        
        // 3. Compute k_2 (spatial interpolation at half step, temporal interp at half step)
        // Estimated half-step location using Euler
        //System.out.println("Half step (k1 -> k2)");
        double[] k2 = stepAhead(
                new double[]{this.getLocation()[0]+k1[0]/2.0,this.getLocation()[1]+k1[1]/2.0},
                elemPart,dep,1.0/2.0,
                neighbours,uvnode,nodexy,trinodes,allelems,u,v,
                tt,st,dt,stepsPerStep,numLayers);

        // 4. Compute k_3 (spatial interpolation at half step, temporal interp at half step)
        //System.out.println("Half step (k2 -> k3)");
        double[] k3 = stepAhead(
                new double[]{this.getLocation()[0]+k2[0]/2.0,this.getLocation()[1]+k2[1]/2.0},
                elemPart,dep,1.0/2.0,
                neighbours,uvnode,nodexy,trinodes,allelems,u,v,
                tt,st,dt,stepsPerStep,numLayers);              
        
        // 5. Compute k_4 (spatial interpolation at end step)
        //System.out.println("End step (k3 -> k4)");
        double[] k4 = stepAhead(
                new double[]{this.getLocation()[0]+k3[0],this.getLocation()[1]+k3[1]},
                elemPart,dep,1.0,
                neighbours,uvnode,nodexy,trinodes,allelems,u,v,
                tt,st,dt,stepsPerStep,numLayers);
        
        // 6. Add it all together
        if (k1[0] == 0 || k2[0] == 0 || k3[0] == 0 || k4[0] == 0)
        {
//            System.out.printf("RK4 attempt to step out of mesh: Location = [%.6e,%.6e], Element = %d\n",
//                    this.getLocation()[0],this.getLocation()[1],elemPart);
        } 
        else 
        {
            advectStep[0] = (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0;
            advectStep[1] = (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6.0;
        }
        return advectStep;
    }
    /**
     * Calculate the correction steps required for the RK4 algorithm
     * @param xy                Location of new velocity to be used (current location plus spatial step)
     * @param elemPart
     * @param dep
     * @param timeStepAhead   Time step ahead i.e. "2.0" if half-step ahead, "1.0" if full step ahead
     * @param neighbours
     * @param uvnode
     * @param nodexy
     * @param trinodes
     * @param allelems
     * @param u
     * @param v
     * @param u1
     * @param v1
     * @param tt
     * @param st
     * @param dt
     * @param stepsPerStep
     * @param fnum
     * @param lastday
     * @param numLayers
     * @return 
     * 
     * 14/02/17 --- REMOVED double[][] u1, double[][] v1, from arguments
     */
    public static double[] stepAhead(double[] xy, int elemPart, int dep, double timeStepAhead,
            int[][] neighbours, double[][] uvnode, double[][] nodexy, int[][] trinodes, int[] allelems,
            double[][] u, double[][] v, 
            int tt, int st, double dt,
            int stepsPerStep, int numLayers)
    {   
        double[] xy_step = new double[2];
        // Generate a "nearest list" for this location
        double[][] xNrList = neighbourCellsList(xy,elemPart,neighbours,uvnode,nodexy,trinodes,allelems);
        // If predicted location of particle is outside mesh, return zero velocity
        if (xNrList[0][0] == 0)
        {
            return xy_step;
        }
        // compute velocities at start and end of entire step, at the new location
        double[] vel = velocityFromNearestList(xNrList,tt,u,v,numLayers,dep);
        // 14/02/17 --- This line changed to directly get velocity from new 7 line velocity files (makes getNextVel redundant)
        double[] velplus1 = velocityFromNearestList(xNrList,tt+1,u,v,numLayers,dep);
        //double[] velplus1 = getNextVel(tt,recordsPerFile,fnum,lastday,xNrList,u,v,u1,v1,numLayers,dep);               
        xy_step[0] = dt*(vel[0] + ((double)(st+1.0/timeStepAhead)/(double)stepsPerStep)*(velplus1[0]-vel[0]));
        xy_step[1] = dt*(vel[1] + ((double)(st+1.0/timeStepAhead)/(double)stepsPerStep)*(velplus1[1]-vel[1]));
        return xy_step;
    }
    
    /**
     * THIS IS THE OLD METHOD FOR 6 ROW VELOCITY FILES
     * 
     * Get the relevant next velocity to compute a time interpolation (which depends on 
     * current time in relation to the records stored in the input arrays) 
     * in the spatial interpolation case
     * @param tt
     * @param recordsPerFile
     * @param fnum
     * @param lastday
     * @param nrList
     * @param u
     * @param v
     * @param u1
     * @param v1
     * @param numLayers
     * @param depLayer
     * @return 
     */
    public static double[] getNextVel(int tt, int recordsPerFile, int fnum, int lastday, double[][] nrList,
            double[][] u, double[][] v, double[][] u1, double[][] v1, int numLayers, int depLayer)
    {
        double[] velplus1 = new double[2];
        if (tt < recordsPerFile-1)
        {
            velplus1 = velocityFromNearestList(nrList,tt+1,u,v,numLayers,depLayer);
        }
        else
        {
            if (fnum < lastday)
            {
                velplus1 = velocityFromNearestList(nrList,0,u1,v1,numLayers,depLayer);
            }
            // Very last timestep: just fix velocity to be same over entire last hour 
            // (as there is no extra file to read in)
            else
            {
                velplus1 = velocityFromNearestList(nrList,tt,u,v,numLayers,depLayer);
            }
        }
        return velplus1;
    }

    public double[] eulerStepOld(double[][] u, double[][] v, double[][] u1, double[][] v1, // velocities
            int[][] neighbours, double[][] uvnode,                                               // other mesh info
            int tt, int st, double dt,                                          // locate particle in space and time
            int stepsPerStep, int recordsPerFile, int fnum, int lastday, int depthLayers,   // info on simulation length
            boolean spatialInterpolate, boolean timeInterpolate)                            // interpolate or not?
    {
        int elemPart = this.getElem();
        int dep = this.getDepthLayer();
        double[] advectStep = new double[2];
        // Set velocity array here first and fill with the right values. Need an array with velocity "now" 
        // and one with velocity at "next" time step, to be populated before go to time interpolation bit 
        // (which should happen last, after any possible spatial interpolation).
        double water_U = 0;
        double water_V = 0;
        // Moved old methods for spatial/temporal interpolation
        if (spatialInterpolate==false)
        {
            // output from matlab is in rows for each time
            // and then columns for each depth within each element
            water_U = u[tt][elemPart*depthLayers+dep];
            water_V = v[tt][elemPart*depthLayers+dep];
            //System.out.printf("Vel(element %d) = %.4f %.4f\n",elemPart,water_U,water_V);
            if (timeInterpolate == true)
            {
                if (tt < recordsPerFile-1)
                {
                    water_U = u[tt][elemPart*depthLayers+dep] + ((double)st/(double)stepsPerStep)*(u[tt+1][elemPart*10+dep]-u[tt][elemPart*depthLayers+dep]);
                    water_V = v[tt][elemPart*depthLayers+dep] + ((double)st/(double)stepsPerStep)*(v[tt+1][elemPart*10+dep]-v[tt][elemPart*depthLayers+dep]);
                } 
                else
                {
                    //water_U = u[tt][elemPart*depthLayers+dep];
                    //water_V = v[tt][elemPart*depthLayers+dep];
                    if (fnum < lastday)
                    {
                        water_U = u[tt][elemPart*depthLayers+dep] + ((double)st/(double)stepsPerStep)*(u1[0][elemPart*10+dep]-u[tt][elemPart*depthLayers+dep] );
                        water_V = v[tt][elemPart*depthLayers+dep] + ((double)st/(double)stepsPerStep)*(v1[0][elemPart*10+dep]-v[tt][elemPart*depthLayers+dep]);                             
                    }
                    // Very last timestep: just fix velocity to be same over entire last hour
                    else
                    {
                        water_U = u[tt][elemPart*depthLayers+dep];
                        water_V = v[tt][elemPart*depthLayers+dep];
                    }
                }
            } 
        }
        else if (spatialInterpolate == true)
        {
            // Print element values for comparison
            //System.out.printf("Vel(element %d) = %.4f %.4f\n",elemPart,u[tt][elemPart*10+dep],v[tt][elemPart*10+dep]);
//                                int nearest = Particle.nearestCentroid(particles[i].getLocation()[0], particles[i].getLocation()[1], uvnode);
//                                int which = Particle.whichElement(particles[i].getLocation()[0], particles[i].getLocation()[1], allelems, nodexy, trinodes);
            this.setNrListToNeighbourCells(neighbours,uvnode);
//            if (nearest != which)
//            {
//                System.out.printf("Particle "+i+" Elem "+particles[i].getElem()+" nearest "+nearest+" whichElem "+which+"\n");
//                double d1 = Particle.distanceEuclid(particles[i].getLocation()[0],particles[i].getLocation()[1],uvnode[which][0],uvnode[which][1]);
//                double d2 = Particle.distanceEuclid(particles[i].getLocation()[0],particles[i].getLocation()[1],uvnode[nearest][0],uvnode[nearest][1]);
//                System.out.printf("dist(elemPart) = "+d1+" dist(nearest) = "+d2+"\n");
//                double[][] a = particles[i].getNrList();
//                System.out.printf("Neighbour cells: \n"+a[0][0]+" "+a[0][1]+"\n"+a[1][0]+" "+a[1][1]+"\n"+a[2][0]+" "+a[2][1]+"\n"+a[3][0]+" "+a[3][1]+"\n");
//            }
            //particles[i].setNrList(uvnode);
            double[] vel = this.velPart(tt,u,v);
            water_U = vel[0];
            water_V = vel[1];
            if (timeInterpolate == true)
            {
                double[] velplus1 = new double[2];
                if (tt < recordsPerFile-1)
                {
                    velplus1 = this.velPart(tt+1,u,v);
                }
                else
                {
                    //velplus1 = particles[i].velocityFromNearestList(0,u1,v1);
                    if (fnum < lastday)
                    {
                        velplus1 = this.velPart(0,u1,v1);
                    }
                    // Very last timestep: just fix velocity to be same over entire last hour
                    else
                    {
                        velplus1 = this.velPart(tt,u,v);
                    }
                }
                water_U = vel[0] + ((double)st/(double)stepsPerStep)*(velplus1[0]-vel[0]);
                water_V = vel[1] + ((double)st/(double)stepsPerStep)*(velplus1[1]-vel[1]);
            }
        }
//        if (st == 1)
//        {
//            System.out.printf("Vel (particle %d element %d) = %.4f %.4f\n",this.ID,elemPart,water_U,water_V);
//        }
        advectStep[0]=dt*water_U;
        advectStep[1]=dt*water_V;
        
        return advectStep;
    }
    
}