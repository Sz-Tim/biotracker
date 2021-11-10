/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package extUtils;

import particle_track.IOUtils;

/**
 * https://www.sanfoundry.com/java-program-check-whether-given-point-lies-given-polygon/
 * @author SA01TA
 */
public class PointInsidePolygon {
    //This is a java program to check whether a point lies in a polygon or not

    private static class Point
    {
        float x, y;
        Point()
        {}

        Point(float p, float q)
        {
            x = p;
            y = q;
        }
        
        @Override
        public String toString()
        {
            return(this.x+" "+this.y);
        }
        
    }

    private static boolean onSegment(Point p, Point q, Point r)
    {
        return q.x <= Math.max(p.x, r.x) && q.x >= Math.min(p.x, r.x)
                && q.y <= Math.max(p.y, r.y) && q.y >= Math.min(p.y, r.y);
    }

    private static int orientation(Point p, Point q, Point r)
    {
        float val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);

        if (val == 0)
            return 0;
        
        return (val > 0) ? 1 : 2;

    }

    private static boolean doIntersect(Point p1, Point q1, Point p2, Point q2)
    {
        int o1 = orientation(p1, q1, p2);
        int o2 = orientation(p1, q1, q2);
        int o3 = orientation(p2, q2, p1);
        int o4 = orientation(p2, q2, q1);

        if (o1 != o2 && o3 != o4)
            return true;

        if (o1 == 0 && onSegment(p1, p2, q1))
            return true;

        if (o2 == 0 && onSegment(p1, q2, q1))
            return true;

        if (o3 == 0 && onSegment(p2, p1, q2))
            return true;

        return o4 == 0 && onSegment(p2, q1, q2);

    }

    public static boolean isInside(Point[] polygon, int n, Point p)
    {
        int INF = 10000;
        if (n < 3)
            return false;
 
        Point extreme = new Point(INF, p.y);

        int count = 0, i = 0;
        do
        {
            int next = (i + 1) % n;
            if (doIntersect(polygon[i], polygon[next], p, extreme))
            {
                if (orientation(polygon[i], p, polygon[next]) == 0)
                    return onSegment(polygon[i], p, polygon[next]);
 
                count++;
            }
            i = next;
        } while (i != 0);

        return (count & 1) == 1;
    }

 

    public static void main(String[] args)
    {
        Point[] polygon1 = { new Point(0, 0), new Point(10, 0),
                new Point(10, 10), new Point(0, 10) };

        int n = 4;

        Point p = new Point(20, 20);
        System.out.println("Point P(" + p.x + ", " + p.y
                + ") lies inside polygon1: " + isInside(polygon1, n, p));

        p = new Point(5, 5);
        System.out.println("Point P(" + p.x + ", " + p.y
                + ") lies inside polygon1: " + isInside(polygon1, n, p));

        Point[] polygon2 = { new Point(0, 0), new Point(5, 5), new Point(5, 0) };
        n = 3;

        p = new Point(3, 3);
        System.out.println("Point P(" + p.x + ", " + p.y

                + ") lies inside polygon2: " + isInside(polygon2, n, p));

        p = new Point(5, 1);
        System.out.println("Point P(" + p.x + ", " + p.y
                + ") lies inside polygon2: " + isInside(polygon2, n, p));

        p = new Point(8, 1);
        System.out.println("Point P(" + p.x + ", " + p.y
                + ") lies inside polygon2: " + isInside(polygon2, n, p));

        Point[] polygon3 = { new Point(0, 0), new Point(10, 0),
                new Point(10, 10), new Point(0, 10), new Point(5, 5) };

        n = 5;

        p = new Point(-1, 10);
        System.out.println("Point P(" + p.x + ", " + p.y
                + ") lies inside polygon3: " + isInside(polygon3, n, p));
        
        
        // Make a check using the actual FVCOM mesh boundary
        float[][] nodexy = IOUtils.readNetcdfFloat2D("C:\\Users\\SA01TA\\Documents\\particle_track\\WestCOMS_mesh.nc","nodexy",null,null);
        int[] bnode = IOUtils.readNetcdfInteger1D("C:\\Users\\SA01TA\\Documents\\particle_track\\WestCOMS_mesh.nc","boundaryNodesAll");
        
        System.out.println("Length of bnodes: "+bnode.length);
        
        Point[] meshOutline = new Point[bnode.length];
        for (int i = 0; i < bnode.length; i++)
        {
            meshOutline[i] = new Point(nodexy[0][bnode[i]-1],nodexy[1][bnode[i]-1]);
        }
        
        System.out.println("Point on mesh outline: "+meshOutline[0].toString());
        
        Point test1 = new Point(-7,56);
        Point test2 = new Point(-8,56);
        Point test3 = new Point((float)-5.6,(float)55.5);
        Point test4 = new Point((float)-5.25,(float)55.5);
        
        System.out.println("Point P(" + test1.x + ", " + test1.y
                + ") lies inside meshOutline: " + isInside(meshOutline, bnode.length, test1));
        System.out.println("Point P(" + test2.x + ", " + test2.y
                + ") lies inside meshOutline: " + isInside(meshOutline, bnode.length, test2));
        System.out.println("Point P(" + test3.x + ", " + test3.y
                + ") lies inside meshOutline: " + isInside(meshOutline, bnode.length, test3));
        System.out.println("Point P(" + test4.x + ", " + test4.y
                + ") lies inside meshOutline: " + isInside(meshOutline, bnode.length, test4));

    }

}

