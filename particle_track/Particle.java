/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package particle_track;

/**
 *
 * @author tomdude
 */
public class Particle {
    
    // horizontal position
    private double[] xy = new double[2];
    private int elem;
    private double[][] nrList = new double[5][2];
    private double[][] cornerList = new double[3][2];
    // vertical position
    private double z;
    private int depLayer;
    // settlement details
    private double age;
    private double density;
    private double mortRate;
    private boolean arrived;
    private boolean viable;
    
    // create a new particle at a defined location, at the water surface
    public Particle(double xstart, double ystart)
    {
        this.xy[0] = xstart;
        this.xy[1] = ystart;
        this.z = 0;
        this.depLayer = 0;
        this.age = 0;
        this.arrived = false;
        this.viable = false;
        this.mortRate = 0.01; // default hourly rate, based on results of Stein et al. 2005 (copepodid, nauplii rate marginally lower = 0.0078)
        this.density = 1;
    }
    
    public double[] getLocation()
    {
        return this.xy;
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
        System.out.printf("setting depth layer: %d (%f, particle depth = %f)\n",depNearest,localDepth*layers[depNearest],this.z);
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
    public boolean getViable()
    {
        return this.viable;
    }
    public boolean getArrived()
    {
        return this.arrived;
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
        // sets the nearestList (NrList) to be identical to the list of elements which neighbour the element that the particle is actually in
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

    
    public double[] velocityFromNearestList(int tt,double u[][],double v[][])
    {
        double[] velocity = new double[2];
        double[] weights = new double[this.nrList.length];
        double usum=0,vsum=0,sum=0;
        for (int i = 0; i < this.nrList.length; i++)
        {
            if (this.nrList[i][1]!=0)
            {
                weights[i]=1.0/(this.nrList[i][1]*this.nrList[i][1]);
            }
            else
            {
                weights[i]=1;
            }
            
            boolean reducedFile = true;
            int numLayers = 10;
            if (reducedFile==true)
            {
                numLayers = 2;
            }
            //System.out.printf("tt %d i %d",tt,i);
            //System.out.printf(" --- elem %d dist %.4f weight %.4f --- vel = %.4f %.4f\n",(int)this.nrList[i][0],this.nrList[i][1],weights[i],u[tt][(int)this.nrList[i][0]*numLayers+this.depLayer],v[tt][(int)this.nrList[i][0]*numLayers+this.depLayer]);
            usum=usum+weights[i]*u[tt][(int)this.nrList[i][0]*numLayers+this.depLayer];
            vsum=vsum+weights[i]*v[tt][(int)this.nrList[i][0]*numLayers+this.depLayer];
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
                //fprintf('particle in %d\n', element);
                break;
            }

        }
        return which;
    }
    
    public static double distanceEuclid(double x1, double y1, double x2, double y2)
    {
        double dist = Math.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
        return dist;
    }
    
}
