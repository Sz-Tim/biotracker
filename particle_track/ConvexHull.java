/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package particle_track;

/**
 * Code to compute the convex hull of a set of points.
 * <p>
 * Modified from:
 * https://rosettacode.org/wiki/Convex_hull#Java
 * Content is available under GNU Free Documentation License 1.2
 * <p>
 * The convex hull encircles the set of points in an anticlockwise direction.
 *
 * @author SA01TA
 */

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
//import java.awt.Point;

import static java.util.Collections.emptyList;

public class ConvexHull {
    private static class MyPoint {
        private final float x;
        private final float y;

        public MyPoint(float x, float y) {
            this.x = x;
            this.y = y;
        }

        //@Override
        public int compareTo(MyPoint o) {
            return Float.compare(x, o.x);
        }

        @Override
        public String toString() {
            return String.format("(%f, %f)", x, y);
        }

        public float[] toFloat() {
            float[] xy = new float[]{this.x, this.y};
            return xy;
        }
    }

    public static float[][] convexHull(float[][] p) {
        // Convert array to a list of points
        List<MyPoint> pt = new ArrayList<>();
        for (float[] floats : p) {
            pt.add(new MyPoint(floats[0], floats[1]));
        }
        // Find the convex hull
        List<MyPoint> hull = convexHull(pt);
        // Convert the list of points back to an array
        float[][] h = new float[hull.size()][2];
        for (int i = 0; i < hull.size(); i++) {
            h[i] = hull.get(i).toFloat();
        }
        return h;
    }

    private static List<MyPoint> convexHull(List<MyPoint> p) {
        if (p.isEmpty()) return emptyList();
        p.sort(MyPoint::compareTo);
        List<MyPoint> h = new ArrayList<>();

        // lower hull
        for (MyPoint pt : p) {
            while (h.size() >= 2 && !ccw(h.get(h.size() - 2), h.get(h.size() - 1), pt)) {
                h.remove(h.size() - 1);
            }
            h.add(pt);
        }

        // upper hull
        int t = h.size() + 1;
        for (int i = p.size() - 1; i >= 0; i--) {
            MyPoint pt = p.get(i);
            while (h.size() >= t && !ccw(h.get(h.size() - 2), h.get(h.size() - 1), pt)) {
                h.remove(h.size() - 1);
            }
            h.add(pt);
        }

        h.remove(h.size() - 1);
        return h;
    }

    // ccw returns true if the three points make a counter-clockwise turn
    private static boolean ccw(MyPoint a, MyPoint b, MyPoint c) {
        return ((b.x - a.x) * (c.y - a.y)) > ((b.y - a.y) * (c.x - a.x));
    }

    public static void main(String[] args) {
        List<MyPoint> points = Arrays.asList(new MyPoint(16, 3),
                new MyPoint(12, 17),
                new MyPoint(0, 6),
                new MyPoint(-4, -6),
                new MyPoint(16, 6),

                new MyPoint(16, -7),
                new MyPoint(16, -3),
                new MyPoint(17, -4),
                new MyPoint(5, 19),
                new MyPoint(19, -8),

                new MyPoint(3, 16),
                new MyPoint(12, 13),
                new MyPoint(3, -4),
                new MyPoint(17, 5),
                new MyPoint(-3, 15),

                new MyPoint(-3, -9),
                new MyPoint(0, 11),
                new MyPoint(-9, -3),
                new MyPoint(-4, -2),
                new MyPoint(12, 10));

        List<MyPoint> hull = convexHull(points);
        System.out.printf("Convex Hull: %s\n", hull);

        float[][] p = new float[points.size()][2];
        for (int i = 0; i < points.size(); i++) {
            p[i] = points.get(i).toFloat();
        }

        float[][] h = convexHull(p);
        for (float[] floats : h) {
            System.out.println(floats[0] + " " + floats[1]);
        }

    }
}