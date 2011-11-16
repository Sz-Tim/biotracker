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
    
    private double[] xy = new double[2];
    private int elem;
    private boolean arrived;
    private boolean viable;
    
    public Particle(double xstart, double ystart)
    {
        this.xy[0] = xstart;
        this.xy[1] = ystart;
        this.arrived = false;
        this.viable = false;
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
    
    public static int nearestCentroid(double x, double y, double[][] uvnode)
    {
        int nearest = 0;
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
        return nearest;
    }
    
    public static int whichElement(double x, double y, int[] elems, double[][] nodexy, int[][] trinodes)
    {
        int which = 0;
        int res = 0;
        for (int i = 0; i < elems.length; i++)
        {
            double[] xt = new double[3];
            double[] yt = new double[3];
            
            for (int j = 0; j < 3; j++)
            {
                 xt[j]=nodexy[trinodes[elems[i]][j]-1][0];
                 yt[j]=nodexy[trinodes[elems[i]][j]-1][1];
            }
            // check whether (x,y) lies within this
            //fprintf('check %d\n', possibleElems(i));
            
            double f1 = (y-yt[0])*(xt[1]-xt[0]) - (x-xt[0])*(yt[1]-yt[0]);
            double f2 = (y-yt[2])*(xt[0]-xt[2]) - (x-xt[2])*(yt[0]-yt[2]);
            double f3 = (y-yt[2])*(xt[2]-xt[1]) - (x-xt[1])*(yt[2]-yt[1]);
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
    
    
    public static int setParticleDepth(int behaviour, int time)
    {
        return 0;
    }
    
}
