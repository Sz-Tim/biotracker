/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package particle_track;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.LineNumberReader;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;

import java.util.*;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.filefilter.WildcardFileFilter;
import ucar.ma2.*;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;

import extUtils.*;

// import netcdf reading libraries for alternative data import
// https://publicwiki.deltares.nl/display/FEWSDOC/Reading+and+writing+netcdf+gridded+data+with+Java
//import ucar.nc2.*;

/**
 * @author sa01ta
 */
public class IOUtils {

    /**
     * Count the lines in a file, in order to read files of unknown length
     *
     * @param aFile
     * @return
     * @throws IOException
     */
    public static int countLines(File aFile) throws Exception {
        LineNumberReader reader = null;
        try {
            reader = new LineNumberReader(new FileReader(aFile));
            while ((reader.readLine()) != null) ;
            return reader.getLineNumber();
        } catch (Exception ex) {
            return -1;
        } finally {
            if (reader != null)
                reader.close();
        }
    }

    public static int countWords(String s) {

        String trim = s.trim();
        if (trim.isEmpty()) {
            return 0;
        }
        return trim.split("\\w+").length;
    }

//    public static double[][] setupStartLocs(String filename, String sitedir, boolean makeLocalCopy)
//    {
//        double startlocs[][] = new double[10][3];
//        System.out.println("Startlocations defined in file: "+filename);
//        File file = new File(sitedir+filename);
//        int nLines = 0;
//        try
//        {
//            nLines = countLines(file);
//            System.out.println("FILE "+file+" NLINES "+nLines);
//        }
//        catch (Exception e)
//        {
//            System.out.println("Cannot open "+sitedir+filename);
//        }
//        startlocs = IOUtils.readFileDoubleArray(sitedir+filename,nLines,5," ",true);
//        
//        if (makeLocalCopy == true)
//        {
//            IOUtils.writeFloatArrayToFile((startlocs, "startlocs.dat",true); 
//        }
//        
//        return startlocs;
//    }

    public static List<HabitatSite> createHabitatSites(String filename, String sitedir, int scaleCol, boolean makeLocalCopy, List<Mesh> meshes, RunProperties rp) {
        List<HabitatSite> habitat = new ArrayList<>();
        System.out.println("Habitat defined in file: " + filename);
        File file;

        if (sitedir != null) {
            file = new File(sitedir + filename);
        } else {
            file = new File(filename);
        }

        int nLines = 0;
        try {
            nLines = countLines(file);
            System.out.println("FILE " + file + " NLINES " + nLines);
        } catch (Exception e) {
            System.err.println("Cannot open " + sitedir + filename);
        }

        // try to create habitat from file
        try {
            BufferedReader in = new BufferedReader(new FileReader(file));    //reading files in specified directory

            String line;
            //boolean printWarning=true;
            int count = 0, countNotCreated = 0;

            while ((line = in.readLine()) != null)    //file reading
            {
                //int numEntries = countWords(line);
                //System.out.println("Creating habitat site "+count);
                String[] values = line.split("\t");
                String ID = values[0];
                float x = (float) Double.parseDouble(values[1]);
                float y = (float) Double.parseDouble(values[2]);
                float depth = 0;
                if (values.length > 3) {
                    depth = (float) Double.parseDouble(values[3]);
                }
                float scale = 0;
                if (values.length > 4 && scaleCol < values.length) {
                    scale = (float) Double.parseDouble(values[scaleCol]);
                }
                //System.out.println("Read locations, need to add site");
                HabitatSite site = new HabitatSite(ID, x, y, depth, scale, meshes, rp);
                //System.out.println("meshType = "+site.getContainingMeshType());
                if (!site.getContainingMeshType().equalsIgnoreCase("NONE")) {
                    habitat.add(site);
                    count++;
                } else {
                    countNotCreated++;
                }

            }
            in.close();
            System.out.println("Created " + count + " habitat sites, skipped " + countNotCreated);
        } catch (Exception e) {
            System.err.println("Cannot create habitat sites from " + file);
        }

        // Make a copy if required
        if (makeLocalCopy) {
            String outName = "startlocs.dat";
            try {
                Files.copy(file.toPath(), Paths.get(outName));
            } catch (Exception e) {
                System.err.println("Cannot copy " + file + " to " + outName);
            }
        }

        return habitat;
    }


    /**
     * Add some extra locations at which settlement is possible, or limit settlement to
     * a smaller selection
     *
     * @param habitat
     * @param sitedir
     * @param startlocs
     * @param limit     // the maximum index of "startlocs" to allow as an endlocation
     * @return
     */
    public static double[][] setupEndLocs(String habitat, String sitedir, double[][] startlocs, int limit) {
        double[][] endlocs = new double[10][3];

        if (habitat.equalsIgnoreCase("fishfarm3_new")) {
            double[][] extraEndLocs = IOUtils.readFileDoubleArray(sitedir + "160119_fishfarms_jun15.dat", 205, 5, " ", true);
            endlocs = new double[startlocs.length + extraEndLocs.length][];
            for (int i = 0; i < startlocs.length; i++) {
                endlocs[i] = startlocs[i].clone();
            }
            for (int i = startlocs.length; i < startlocs.length + endlocs.length; i++) {
                endlocs[i] = extraEndLocs[i - startlocs.length].clone();
            }
        } else if (limit == 0) {
            endlocs = new double[startlocs.length][];
            for (int i = 0; i < startlocs.length; i++) {
                endlocs[i] = startlocs[i].clone();
            }
        } else {
            endlocs = new double[limit][];
            for (int i = 0; i < limit; i++) {
                endlocs[i] = startlocs[i].clone();
            }
        }
        System.out.println("Endlocs NLINES " + endlocs.length);
        return endlocs;
    }

//    public static double[][] setupOpenBCLocs(String location, String datadir2)
//    {
//        double open_BC_locs[][] = new double[10][3];
//        if (location.equalsIgnoreCase("minch") || location.equalsIgnoreCase("minch_continuous") || location.equalsIgnoreCase("minch_jelly"))
//        {
//            //String sitedir=basedir+"minch_sites/";
//            //open_BC_locs = IOUtils.readFileDoubleArray(sitedir+"open_boundary_locs.dat",120,3," ",true);
//            open_BC_locs = IOUtils.readFileDoubleArray(datadir2+"open_boundary_locs_os.dat",137,3," ",true);
//        } 
//        else
//        {
//            open_BC_locs = IOUtils.readFileDoubleArray("C:\\Users\\sa01ta\\Documents\\lorn\\120903_renewableimpact\\131107_revision\\open_boundary_locs.dat",93,3," ",true);
//        }
//        return open_BC_locs;
//    }

    public static List<Particle> readRestartParticles(RunProperties rp) {
        List<Particle> parts = new ArrayList<Particle>();
        System.out.println("Particles to restart defined in file: " + rp.restartParticles);
        File file = new File(rp.restartParticles);
        int nLines = 0;
        try {
            nLines = countLines(file);
            System.out.println("FILE " + file + " NLINES " + nLines);
        } catch (Exception e) {
            System.err.println("Cannot open " + rp.restartParticles);
        }

        int[] todayInt = ISO_datestr.dateIntParse(rp.start_ymd);
        ISO_datestr today = new ISO_datestr(todayInt[0], todayInt[1], todayInt[2]);
        int todayDateNum = today.getDateNum();
        //System.out.println("----- In readRestartParticles, date = "+today.getDateStr()+" ("+todayDateNum+") -----");

        try {
            BufferedReader in = new BufferedReader(new FileReader(file));    //reading files in specified directory

            String line;
            //boolean printWarning=true;
            int count = 0;

            while ((line = in.readLine()) != null)    //file reading
            {
                // Ignore the header line
                if (count > 0) {
                    //int numEntries = countWords(line);
                    //System.out.println("Creating restart particle "+count);
                    Particle p = new Particle(line, rp.species);
                    boolean add = true;
                    if (rp.restartParticlesCutoffDays > 0) {
                        ISO_datestr pStart = p.getStartDate();
                        int pStartDateNum = pStart.getDateNum();
                        //System.out.println("particle startDate = "+pStart.getDateStr()+" ("+pStartDateNum+")");
                        if (todayDateNum - pStartDateNum > rp.restartParticlesCutoffDays) {
                            //System.out.println("Particle NOT added");
                            add = false;
                        }
                    }
                    if (add) {
                        //System.out.println("Particle added");
                        parts.add(p);
                    }
                }
                count++;
            }
            in.close();
            System.out.println("Created " + parts.size() + " particles");
        } catch (Exception e) {
            System.err.println("Cannot create particles from " + file);
        }

        return parts;
    }

    public static int[] readFileInt1D(String fullFileName) throws Exception {
        //double open_BC_locs[][] = new double[10][3];
        int nLines = countLines(new File(fullFileName));
        int[][] vals = IOUtils.readFileIntArray(fullFileName, nLines, 1, " ", true);
        int[] valsOut = new int[nLines];
        for (int i = 0; i < nLines; i++) {
            valsOut[i] = vals[i][0];
        }
        return valsOut;
    }


    public static double[] readFileDouble1D(String fullFileName) throws Exception {
        int nLines = countLines(new File(fullFileName));
        double[][] vals = IOUtils.readFileDoubleArray(fullFileName, nLines, 1, " ", false);
        double[] valsOut = new double[nLines];
        for (int i = 0; i < nLines; i++) {
            valsOut[i] = vals[i][0];
        }
        return valsOut;
    }


    public static float[] readNetcdfFloat1D(String filename, String variable) {
        System.out.println("Reading variable: " + variable);
        float[] floatOut = null;

        try (NetcdfFile dataFile = NetcdfFile.open(filename, null)) {
            Variable dataVar = dataFile.findVariable(variable);
            if (dataVar == null) {
                System.out.println("  Can't find Variable: " + variable);
                return (floatOut);
            }

            int[] shape = dataVar.getShape();
            //int origin = 0;
            if (shape.length > 1) {
                return (floatOut);
            }

            try {
                ArrayFloat.D1 dataArray = (ArrayFloat.D1) dataVar.read();
                // Put the values into a native array
                floatOut = new float[shape[0]];
                for (int d1 = 0; d1 < shape[0]; d1++) {
                    floatOut[d1] = dataArray.get(d1);
                }
            } catch (ClassCastException e) {
                ArrayDouble.D1 dataArray = (ArrayDouble.D1) dataVar.read();
                // Put the values into a native array
                floatOut = new float[shape[0]];
                for (int d1 = 0; d1 < shape[0]; d1++) {
                    floatOut[d1] = (float) dataArray.get(d1);
                }
            }


        } catch (IOException ioe) {
            ioe.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        } finally {

        }
        return floatOut;
    }

    public static int[] readNetcdfInteger1D(String filename, String variable) {
        System.out.println("Reading variable: " + variable);
        int[] intOut = null;

        try (NetcdfFile dataFile = NetcdfFile.open(filename, null)) {
            Variable dataVar = dataFile.findVariable(variable);
            if (dataVar == null) {
                System.out.println("  Can't find Variable: " + variable);
                return (intOut);
            }

            int[] shape = dataVar.getShape();
            //int origin = 0;
            ArrayInt.D1 dataArray = (ArrayInt.D1) dataVar.read();
            // Put the values into a native array
            intOut = new int[shape[0]];
            for (int d1 = 0; d1 < shape[0]; d1++) {
                intOut[d1] = dataArray.get(d1);
            }
        } catch (IOException ioe) {
            ioe.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return intOut;
    }

    /**
     * Method to read NetCDF variable
     *
     * @param filename
     * @param variable
     * @param origin
     * @param shape
     * @return
     * @throws IOException
     * @throws InvalidRangeException
     */
    public static float[][] readNetcdfFloat2D(String filename, String variable, int[] origin, int[] shape) {
//        System.out.println("Reading variable: "+variable);
        //float[][] floatOut = new float[2][10];
        float[][] floatOut = null;

        try (NetcdfFile dataFile = NetcdfFile.open(filename, null)) {
            Variable dataVar = dataFile.findVariable(variable);
            if (dataVar == null) {
                System.err.println("  Can't find Variable: " + variable);
                return (floatOut);
            }
//            int[] origin = new int[2];
//            int[] shape = dataVar.getShape();

            // Check that origin is within the array
            if (origin == null
                    || origin[0] > dataVar.getShape()[0]
                    || origin[1] > dataVar.getShape()[1]) {
//                System.out.println("Origin not supplied, or outside bounds of variable "+variable);
                origin = new int[3];
            }

            //int[] shape = dataVar.getShape();
            // Check that origin+shape is within the array
            if (shape == null
                    || origin[0] + shape[0] > dataVar.getShape()[0]
                    || origin[1] + shape[1] > dataVar.getShape()[1]) {
//                System.out.println("Shape not supplied, or outside bounds of variable "+variable);
                shape = dataVar.getShape();
            }

            System.out.println(variable + " (" + shape[0] + "," + shape[1] + ")");

            try {
                ArrayFloat.D2 dataArray = (ArrayFloat.D2) dataVar.read(origin, shape);
                // Put the values into a native array
                floatOut = new float[shape[0]][shape[1]];
                for (int d1 = 0; d1 < shape[0]; d1++) {
                    for (int d2 = 0; d2 < shape[1]; d2++) {
                        floatOut[d1][d2] = dataArray.get(d1, d2);
                    }
                }
            } catch (ClassCastException e) {
                ArrayDouble.D2 dataArray = (ArrayDouble.D2) dataVar.read(origin, shape);
                // Put the values into a native array
                floatOut = new float[shape[0]][shape[1]];
                for (int d1 = 0; d1 < shape[0]; d1++) {
                    for (int d2 = 0; d2 < shape[1]; d2++) {
                        floatOut[d1][d2] = (float) dataArray.get(d1, d2);
                    }
                }
            }

        } catch (IOException ioe) {
            ioe.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return floatOut;
    }

    public static int[][] readNetcdfInteger2D(String filename, String variable) {
        System.out.println("Reading variable: " + variable);
        //int[][] intOut = new int[2][10];
        int[][] intOut = null;

        try (NetcdfFile dataFile = NetcdfFile.open(filename, null)) {
            Variable dataVar = dataFile.findVariable(variable);
            if (dataVar == null) {
                System.out.println("  Can't find Variable: " + variable);
                return (intOut);
            }

            int[] shape = dataVar.getShape();
            int[] origin = new int[2];
            ArrayInt.D2 dataArray = (ArrayInt.D2) dataVar.read(origin, shape);
            // Put the values into a native array
            intOut = new int[shape[0]][shape[1]];
            for (int d1 = 0; d1 < shape[0]; d1++) {
                for (int d2 = 0; d2 < shape[1]; d2++) {
                    intOut[d1][d2] = dataArray.get(d1, d2);
                }
            }
        } catch (IOException ioe) {
            ioe.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return intOut;
    }

    /**
     * Reshape a two dimensional array of floating point numbers
     *
     * @param A
     * @param m
     * @param n
     * @return
     */
    public static float[][] reshapeFloat(float[][] A, int m, int n) {
        int origM = A.length;
        int origN = A[0].length;
        if (origM * origN != m * n) {
            throw new IllegalArgumentException("New matrix must be of same area as matix A");
        }
        float[][] B = new float[m][n];
        float[] A1D = new float[A.length * A[0].length];

        int index = 0;
        for (float[] floats : A) {
            for (int j = 0; j < A[0].length; j++) {
                A1D[index++] = floats[j];
            }
        }

        index = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                B[j][i] = A1D[index++];
            }

        }
        return B;
    }

    /**
     * @param filename
     * @param variable
     * @param origin
     * @param shape
     * @return
     * @throws IOException
     * @throws InvalidRangeException
     */
    public static float[][][] readNetcdfFloat3D(String filename, String variable, int[] origin, int[] shape) {
//        System.out.println("Reading variable: "+variable);
        //float[][][] floatOut = new float[2][2][10];
        float[][][] floatOut = null;

        try (NetcdfFile dataFile = NetcdfFile.open(filename, null)) {
            Variable dataVar = dataFile.findVariable(variable);
            if (dataVar == null) {
                System.err.println("  Can't find Variable: " + variable);
                return (floatOut);
            }

            // set origin externally, always?
            // set shape externally. based on known size of mesh required. Handle case when it doesn't work here, rather than needing to set here?

            //int[] origin = new int[3];
            // Check that origin is within the array
            if (origin == null
                    || origin[0] > dataVar.getShape()[0]
                    || origin[1] > dataVar.getShape()[1]
                    || origin[2] > dataVar.getShape()[2]) {
//                System.out.println("Origin not supplied, or outside bounds of variable "+variable);
                origin = new int[3];
            }

            //int[] shape = dataVar.getShape();
            // Check that origin+shape is within the array
            if (shape == null
                    || origin[0] + shape[0] > dataVar.getShape()[0]
                    || origin[1] + shape[1] > dataVar.getShape()[1]
                    || origin[2] + shape[2] > dataVar.getShape()[2]) {
//                System.out.println("Shape not supplied, or outside bounds of variable "+variable);
                shape = dataVar.getShape();
            }

            System.out.println(variable + " (" + shape[0] + "," + shape[1] + "," + shape[2] + ")");

            try {
                ArrayFloat.D3 dataArray = (ArrayFloat.D3) dataVar.read(origin, shape);
                // Put the values into a native array
                floatOut = new float[shape[0]][shape[1]][shape[2]];
                for (int d1 = 0; d1 < shape[0]; d1++) {
                    for (int d2 = 0; d2 < shape[1]; d2++) {
                        for (int d3 = 0; d3 < shape[2]; d3++) {
                            floatOut[d1][d2][d3] = dataArray.get(d1, d2, d3);
                        }
                    }
                }
            } catch (ClassCastException e) {
                ArrayDouble.D3 dataArray = (ArrayDouble.D3) dataVar.read(origin, shape);
                // Put the values into a native array
                floatOut = new float[shape[0]][shape[1]][shape[2]];
                for (int d1 = 0; d1 < shape[0]; d1++) {
                    for (int d2 = 0; d2 < shape[1]; d2++) {
                        for (int d3 = 0; d3 < shape[2]; d3++) {
                            floatOut[d1][d2][d3] = (float) dataArray.get(d1, d2, d3);
                        }
                    }
                }
            }

        } catch (IOException ioe) {
            ioe.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return floatOut;
    }

    /**
     * @param filename
     * @param variable
     * @param origin
     * @param shape
     * @return
     * @throws IOException
     * @throws InvalidRangeException
     */
    public static float[][][][] readNetcdfFloat4D(String filename, String variable, int[] origin, int[] shape) {
        System.out.println("Reading variable: " + variable);
        //float[][][][] floatOut = new float[2][2][10][1];
        float[][][][] floatOut = null;

        try (NetcdfFile dataFile = NetcdfFile.open(filename, null)) {
            Variable dataVar = dataFile.findVariable(variable);
            if (dataVar == null) {
                System.err.println("  Can't find Variable: " + variable);
                return (floatOut);
            }

            // set origin externally, always?
            // set shape externally. based on known size of mesh required. Handle case when it doesn't work here, rather than needing to set here?

            //int[] origin = new int[3];
            // Check that origin is within the array
            if (origin == null
                    || origin[0] > dataVar.getShape()[0]
                    || origin[1] > dataVar.getShape()[1]
                    || origin[2] > dataVar.getShape()[2]
                    || origin[3] > dataVar.getShape()[3]) {
                System.out.println("Origin not supplied, or outside bounds of variable " + variable);
                origin = new int[3];
            }

            //int[] shape = dataVar.getShape();
            // Check that origin+shape is within the array
            if (shape == null
                    || origin[0] + shape[0] > dataVar.getShape()[0]
                    || origin[1] + shape[1] > dataVar.getShape()[1]
                    || origin[2] + shape[2] > dataVar.getShape()[2]
                    || origin[3] + shape[3] > dataVar.getShape()[3]) {
                System.out.println("Shape not supplied, or outside bounds of variable " + variable);
                shape = dataVar.getShape();
            }

            System.out.println(variable + " (" + shape[0] + "," + shape[1] + "," + shape[2] + "," + shape[3] + ")");

            try {
                ArrayFloat.D4 dataArray = (ArrayFloat.D4) dataVar.read(origin, shape);
                // Put the values into a native array
                floatOut = new float[shape[0]][shape[1]][shape[2]][shape[3]];
                for (int d1 = 0; d1 < shape[0]; d1++) {
                    for (int d2 = 0; d2 < shape[1]; d2++) {
                        for (int d3 = 0; d3 < shape[2]; d3++) {
                            for (int d4 = 0; d4 < shape[2]; d4++) {
                                floatOut[d1][d2][d3][d4] = dataArray.get(d1, d2, d3, d4);
                            }
                        }
                    }
                }
            } catch (ClassCastException e) {
                ArrayDouble.D4 dataArray = (ArrayDouble.D4) dataVar.read(origin, shape);
                // Put the values into a native array
                floatOut = new float[shape[0]][shape[1]][shape[2]][shape[3]];
                for (int d1 = 0; d1 < shape[0]; d1++) {
                    for (int d2 = 0; d2 < shape[1]; d2++) {
                        for (int d3 = 0; d3 < shape[2]; d3++) {
                            for (int d4 = 0; d4 < shape[2]; d4++) {
                                floatOut[d1][d2][d3][d3] = (float) dataArray.get(d1, d2, d3, d3);
                            }
                        }
                    }
                }
            }

        } catch (IOException ioe) {
            ioe.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return floatOut;
    }


    public static double[][] readFileDoubleArray(String filename, int rows, int cols, String sep, boolean note) {
        double[][] myDouble = new double[rows][cols];
        int x = 0, y = 0;
        boolean failed = false;
        double sum = 0;
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));    //reading files in specified directory

            String line;
            boolean printWarning = true;

            outer:
            while ((line = in.readLine()) != null)    //file reading
            {

                int numEntries = countWords(line);
                if (numEntries < cols && printWarning) {
                    System.out.println("WARNING: Number of entries on line = " + numEntries);
                    printWarning = false;
                }
                y = 0;
                if (x >= rows) {
                    System.out.println(filename + " has more rows than expected.");
                    failed = true;
                    break;
                }
                String[] values = line.split(" ");
                for (String str : values) {
                    if (y >= cols) {
                        System.out.println(filename + " has more columns than expected: " + y);
                        failed = true;
                        break outer;
                    }
                    double str_double = Double.parseDouble(str);
                    myDouble[x][y] = str_double;
                    sum += myDouble[x][y];
                    //System.out.print(myDouble[x][y] + " ");
                    y++;

                }

                //System.out.println("");
                x++;
            }
            in.close();

        } catch (IOException ioException) {
            throw new RuntimeException(ioException);
            //System.err.println("******************* Cannot read from file "+filename+" ******************************");
            //failed = true;
        }
        if (note && !failed) {
            System.out.printf("Created %dx%d array from file: %s\n", myDouble.length, myDouble[0].length, filename);
            System.out.println("Array sum at read time = " + sum);
        } else if (failed) {
            System.out.println("FAILED to read file " + filename);
            System.exit(1);
        }
        return myDouble;
    }

    public static int[][] readFileIntArray(String filename, int rows, int cols, String sep, boolean note) {
        int[][] myInt = new int[rows][cols];
        int x = 0, y = 0;
        boolean failed = false;
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));    //reading files in specified directory
            //System.out.println("in readFileIntArray");
            String line;
            while ((line = in.readLine()) != null)    //file reading
            {
                y = 0;
                if (x >= rows) {
                    System.out.println("WARNING: " + filename + " has more rows than expected. EXITING!");
                    failed = true;
                    break;
                }
                String[] values = line.split(sep);
                for (String str : values) {
                    if (y >= cols) {
                        System.out.println("WARNING: " + filename + " has more columns than expected. Some information not read!");
                        //failed = true;
                        //break;
                    }
                    int str_int = Integer.parseInt(str);
                    myInt[x][y] = str_int;
                    //System.out.print(myDouble[x][y] + " ");
                    y++;

                }
                //System.out.println("");
                x++;
            }
            in.close();


        } catch (Exception e) {
            System.err.println("******************* Cannot read from file " + filename + " ******************************");
            failed = true;
        }
        if (note && !failed) {
            System.out.printf("Created %dx%d array from file: %s\n", myInt.length, myInt[0].length, filename);
        } else if (failed) {
            System.out.println("FAILED to read file " + filename);
            System.exit(1);
        }
        return myInt;
    }

    /**
     * Reads in sunrise and sunset times rounded to the nearest hour from the specified file. File should have the first
     * column as a YYYYMMDD date, second column as the sunrise hour, and the third column as the sunset hour. Dates are
     * trimmed based on the simulation start and end dates.
     *
     * @param filename full path to daylight data file
     * @param startDate ISO_datestr starting date
     * @param endDate ISO_datestr ending date
     * @param numberOfDays number of simulation days
     * @param sep text file delimiter
     * @param note print messages?
     *
     */
    public static int[][] readDaylightHours(String filename, ISO_datestr startDate, ISO_datestr endDate, int numberOfDays, String sep, boolean note) {
        int[][] daylightHours = new int[numberOfDays][2];
        int dayCount = 0;
        boolean failed = false;
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));    //reading files in specified directory

            String line;
            boolean printWarning = true;

            outer:
            while ((line = in.readLine()) != null) {
                int numEntries = countWords(line);
                if (numEntries < 3 && printWarning) {
                    System.out.println("WARNING: Number of entries on line = " + numEntries);
                    printWarning = false;
                }
                String[] values = line.split(" ");
                ISO_datestr rowDate = new ISO_datestr(values[0]);
                if (rowDate.getDateNum() >= startDate.getDateNum() && rowDate.getDateNum() < endDate.getDateNum()) {
                    if (note) {
                        System.out.println("Adding sunrise & sunset for " + rowDate);
                    }
                    daylightHours[dayCount][0] = Integer.parseInt(values[1]);
                    daylightHours[dayCount][1] = Integer.parseInt(values[2]);
                    dayCount++;
                }

            }
            in.close();

        } catch (IOException ioException) {
            throw new RuntimeException(ioException);
            //System.err.println("******************* Cannot read from file "+filename+" ******************************");
            //failed = true;
        }
        if (note && !failed) {
            System.out.printf("Created %dx%d array from file: %s\n", daylightHours.length, daylightHours[0].length, filename);
        } else if (failed) {
            System.out.println("FAILED to read file " + filename);
            System.exit(1);
        }
        return daylightHours;
    }

    public static void writeMovementsHeader(String header, String filename) {
        try {
            FileWriter fstream = new FileWriter(filename, false);
            PrintWriter out = new PrintWriter(fstream);
            out.println(header);
            out.close();
        } catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
        }
    }
    public static void writeMovements(Particle part, int currentHour, int st, double[] displacement, double[] advectStep, double[] activeMovement, double[] diffusion, String filename, boolean append) {
        try {
            // Create file
            FileWriter fstream = new FileWriter(filename, append);
            PrintWriter out = new PrintWriter(fstream);
            out.printf("%d %d %d %s %.1f %.4f %.4f %.4f %d %d %.4f",
                    currentHour,
                    st,
                    part.getID(),
                    part.getStartDate().getDateStr(),
                    part.getAge(),
                    part.getLocation()[0],
                    part.getLocation()[1],
                    part.getDepth(),
                    part.getDepthLayer(),
                    part.getStatus(),
                    part.getDegreeDays()
                    );
            for (double dim: displacement) {
                out.printf(" %.4f", dim);
            }
            for (double dim: advectStep) {
                out.printf(" %.4f", dim);
            }
            for (double dim: activeMovement) {
                out.printf(" %.4f", dim);
            }
            for (double dim: diffusion) {
                out.printf(" %.4f", dim);
            }
            out.print("\n");
            //Close the output stream
            out.close();
        } catch (Exception e) {//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }

    }


    /**
     * Print a predetermined string of characters to a file. Also ensures a new
     * file is started for a given day for particle locations
     *
     * @param headerString
     * @param filename
     */
    public static void printFileHeader(String headerString, String filename) {
        try {
            FileWriter fstream = new FileWriter(filename);
            PrintWriter out = new PrintWriter(fstream);
            out.printf(headerString);
            out.printf("\n");
            out.close();
        } catch (Exception e) {//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }

    /**
     * Print an array of doubles directly to a file,
     *
     * @param variable
     * @param filename
     * @param asInt
     */
    public static void writeFloatArrayToFile(float[][] variable, String filename, boolean asInt, boolean firstColInt) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(filename);
            PrintWriter out = new PrintWriter(fstream);
            for (int i = 0; i < variable.length; i++) {
                for (int j = 0; j < variable[0].length; j++) {
                    if (asInt || (firstColInt && j == 0)) {
                        out.printf("%d", (int) variable[i][j]);
                    } else {
                        out.printf("%.4e", variable[i][j]);
                    }
                    if (j < variable[0].length - 1) {
                        out.printf(" ");
                    }
                }
                out.printf("\n");
            }
            //Close the output stream
            out.close();
        } catch (Exception e) {//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }

    //    /**
//     * As previous method but as floats
//     * @param variable
//     * @param filename 
//     */
//    public static void writeDoubleArrayToFile2(double[][] variable, String filename)
//    {
//        try
//        {
//            // Create file 
//            FileWriter fstream = new FileWriter(filename);
//            PrintWriter out = new PrintWriter(fstream);
//            for (int i = 0; i < variable.length; i++)
//            {
//                for (int j = 0; j < variable[0].length; j++)
//                {
//                    out.printf("%f ",variable[i][j]);
//                }
//                out.printf("\n");
//            }
//            //Close the output stream
//            out.close();
//        }catch (Exception e){//Catch exception if any
//            System.err.println("Error: " + e.getMessage());
//        }
//    }
    public static void writeIntegerArrayToFile(int[][] variable, String filename) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(filename);
            PrintWriter out = new PrintWriter(fstream);
            for (int i = 0; i < variable.length; i++) {
                for (int j = 0; j < variable[0].length; j++) {
                    out.printf("%d ", variable[i][j]);
                }
                out.printf("\n");
            }
            //Close the output stream
            out.close();
        } catch (Exception e) {//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }

    /**
     * Write (append, including time) particle locations to file.
     * One particle for each source site.
     * Allows plotting of single trajectory from each source site.
     *
     * @param particles
     * @param npartsSaved the number of particles for which to save locations (from id=0)
     * @param tt
     * @param filename
     */
    public static void particleLocsToFile(Particle[] particles, int npartsSaved, int tt, String filename) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(filename, true);
            PrintWriter out = new PrintWriter(fstream);
            for (int i = 0; i < npartsSaved; i++) {
                out.printf("%d %d %.1f %.1f %.1f %d %d %d\n", tt, particles[i].getID(),
                        particles[i].getLocation()[0], particles[i].getLocation()[1], particles[i].getDepth(), particles[i].getDepthLayer(),
                        particles[i].getElem(), particles[i].getStatus());
            }
            //Close the output stream
            out.close();
        } catch (Exception e) {//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }

    public static void particleLocsToFile(List<Particle> particles, int npartsSaved, int tt, String filename) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(filename, true);
            PrintWriter out = new PrintWriter(fstream);
            for (int i = 0; i < npartsSaved; i++) {
                out.printf("%d %d %.1f %.1f %.1f %d %d %d\n", tt, particles.get(i).getID(),
                        particles.get(i).getLocation()[0], particles.get(i).getLocation()[1], particles.get(i).getDepth(), particles.get(i).getDepthLayer(),
                        particles.get(i).getElem(), particles.get(i).getStatus());
            }
            //Close the output stream
            out.close();
        } catch (Exception e) {//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }

    /**
     * Print ALL particle locations at a particular time, with corresponding start locations.
     *
     * @param particles
     * @param filename
     */
    public static void particleLocsToFile1(Particle[] particles, String filename) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(filename, true);
            PrintWriter out = new PrintWriter(fstream);
            for (int i = 0; i < particles.length; i++) {
                out.printf("%d %f %f %f %d %f %f\n", i,
                        particles[i].getStartLocation()[0], particles[i].getStartLocation()[1], particles[i].getDepth(), particles[i].getDepthLayer(),
                        particles[i].getLocation()[0], particles[i].getLocation()[1]);
            }
            //Close the output stream
            out.close();
        } catch (Exception e) {//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }

    public static void particleLocsToFile1(List<Particle> particles, String filename, boolean append) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(filename, append);
            PrintWriter out = new PrintWriter(fstream);
            for (int i = 0; i < particles.size(); i++) {
                out.printf("%d %f %f %f %f %f %d\n", i,
                        particles.get(i).getStartLocation()[0], particles.get(i).getStartLocation()[1],
                        particles.get(i).getLocation()[0], particles.get(i).getLocation()[1], particles.get(i).getDepth(), particles.get(i).getDepthLayer());
            }
            //Close the output stream
            out.close();
        } catch (Exception e) {//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }

    /**
     * Write all current particle information to file - new output defined November 2017
     *
     * @param particles
     * @param currentHour
     * @param filename
     * @param append
     */
    public static void particleLocsToFile_full(List<Particle> particles, int currentHour, String filename, boolean append) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(filename, append);
            PrintWriter out = new PrintWriter(fstream);
            for (Particle p : particles) {
                out.printf("%d %d %s %.1f %s %.1f %.1f %.1f %d %d %d %.4f %d\n",
                        currentHour,
                        p.getID(),
                        p.getStartDate().getDateStr(),
                        p.getAge(),
                        p.getStartID(),
                        p.getLocation()[0],
                        p.getLocation()[1],
                        p.getDepth(),
                        p.getDepthLayer(),
                        p.getElem(),
                        p.getStatus(),
                        p.getDensity(),
                        p.getMesh()
                );
            }

            //Close the output stream
            out.close();
        } catch (Exception e) {//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }

    /**
     * @param particles
     * @param currentHour
     * @param filename
     * @param append
     */
    public static void particlesToRestartFile(List<Particle> particles, int currentHour, String filename, boolean append, RunProperties rp) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(filename, append);
            PrintWriter out = new PrintWriter(fstream);
            for (Particle p : particles) {
                if (rp.coordRef.equalsIgnoreCase("WGS84")) {
                    out.printf("%d %d %s %.1f %s %.7f %.7f %d %d %.4f %d %.2f %d %.2f %.2f %.2f %.2f\n",
                            currentHour,
                            p.getID(),
                            p.getStartDate().getDateStr(),
                            p.getAge(),
                            p.getStartID(),
                            p.getLocation()[0],
                            p.getLocation()[1],
                            p.getElem(),
                            p.getStatus(),
                            p.getDensity(),
                            p.getMesh(),
                            p.getDepth(),
                            p.getDepthLayer(),
                            p.getDegreeDays(),
                            p.getxTotal(),
                            p.getyTotal(),
                            p.getzTotal()
                    );
                } else {
                    out.printf("%d %d %s %.1f %s %.1f %.1f %d %d %.4f %d %.2f %d %.2f %.2f %.2f %.2f\n",
                            currentHour,
                            p.getID(),
                            p.getStartDate().getDateStr(),
                            p.getAge(),
                            p.getStartID(),
                            p.getLocation()[0],
                            p.getLocation()[1],
                            p.getElem(),
                            p.getStatus(),
                            p.getDensity(),
                            p.getMesh(),
                            p.getDepth(),
                            p.getDepthLayer(),
                            p.getDegreeDays(),
                            p.getxTotal(),
                            p.getyTotal(),
                            p.getzTotal()
                    );
                }
            }

            //Close the output stream
            out.close();
            System.out.println("\n---- Writing locations for restart file");
        } catch (Exception e) {//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }


    public static void particleLocsToNetcdfFile() {
        // https://www.unidata.ucar.edu/software/netcdf/examples/programs/Simple_xy_wr.java
        // What about appending data? Need to check how that would work.
    }


    public static void arrivalToFile(Particle p, ISO_datestr currentDate, double currentTime, String filename, boolean append) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(filename, append);
            PrintWriter out = new PrintWriter(fstream);
//            for (Particle p : particles)
//            {
            //System.out.printf("Settlement\n");
            out.printf("%d %s %.2f %s %s %.2f %s %f %f\n",
                    p.getID(),
                    p.getStartDate().getDateStr(),
                    p.getStartTime(),
                    p.getStartID(),
                    currentDate.getDateStr(),
                    currentTime,
                    p.getLastArrival(),
                    p.getAge(),
                    p.getDensity()
            );
            //}

            //Close the output stream
            out.close();
        } catch (Exception e) {//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }


    /**
     * Print (append) single particle location to file
     *
     * @param time
     * @param particles
     * @param i
     * @param filename
     */
    public static void particleLocsToFile2(double time, Particle[] particles, int i, String filename) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(filename, true);
            PrintWriter out = new PrintWriter(fstream);
            out.printf("%f %f %f %d\n", time, particles[i].getLocation()[0], particles[i].getLocation()[1], particles[i].getElem());
            //Close the output stream
            out.close();
        } catch (Exception e) {//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }


    /**
     * Make additions to the element presence counts (PSTEPS)
     *
     * @param particles
     * @param rp
     * @param pstepsMature
     * @param pstepsImmature
     * @param dt
     */
    public static void pstepsUpdater(List<Particle> particles, RunProperties rp,
                                     float[][] pstepsMature, float[][] pstepsImmature, double dt) {
        //System.out.println("size pstepsImmature: "+pstepsMature.length+" "+pstepsMature[0].length);
        //System.out.println("size pstepsMature: "+pstepsImmature.length+" "+pstepsImmature[0].length);

        for (Particle p : particles) {
            double d = p.getDensity();
            //System.out.println("density = "+d+" mortRate = "+p.getMortRate());
            int elemPart = p.getElem();

            int col = 0;
            if (rp.splitPsteps) {
                col = p.getStartIndex();
            }

            if (p.getViable()) {
                pstepsMature[elemPart][col] += (float) d * (float) (dt / 3600);//*1.0/rp.stepsPerStep;
            } else if (p.getFree()) {
                //System.out.println("Printing to pstepsImmature");
                pstepsImmature[elemPart][col] += (float) d * (float) (dt / 3600);//*1.0/rp.stepsPerStep;
            }
        }
    }

    public static void pstepsSparseUpdater(List<Particle> particles, RunProperties rp,
                                           List<SparseFloatArray> pstepsMature, List<SparseFloatArray> pstepsImmature, double subStepDt) {
        for (Particle p : particles) {
            double d = p.getDensity();
            int elemPart = p.getElem();
            int col = 1;
            if (rp.splitPsteps) {
                col = p.getStartIndex() + 1;
            }
            if (p.getViable()) {
                SparseFloatArray arr = pstepsMature.get(col); // Get the relevant sparse array from the list
                arr.put(col, arr.get(elemPart) + (float) d * (float) (subStepDt / 3600)); // place the new value in the array
                pstepsMature.set(col, arr); // put the array back in the list

            } else if (p.getFree()) {
                SparseFloatArray arr = pstepsImmature.get(col); // Get the relevant sparse array from the list
                arr.put(col, arr.get(elemPart) + (float) d * (float) (subStepDt / 3600)); // place the new value in the array
                pstepsImmature.set(col, arr); // put the array back in the list
            }
        }
    }

    public static ArrayList<String> checkHydroFilesExist(RunProperties rp, ISO_datestr startDate, ISO_datestr endDate, int numberOfDays) {
        ArrayList<String> missingHydroFiles = new ArrayList<>();
        ISO_datestr checkDate = new ISO_datestr(startDate.getDateStr());
        System.out.println("Checking files from " + startDate.getDateStr() + " to " + endDate.getDateStr());
        for (int i = 0; i < numberOfDays; i++) {
            List<File> checkFile = (List<File>) FileUtils.listFiles(
                    new File(rp.datadir + rp.datadirPrefix + checkDate.getYear() + rp.datadirSuffix + System.getProperty("file.separator")),
                    new WildcardFileFilter(rp.location + rp.minchVersion + "_" + checkDate.getYear() + String.format("%02d", checkDate.getMonth()) + String.format("%02d", checkDate.getDay()) + "*.nc"),
                    null);
            if (checkFile.isEmpty()) {
                missingHydroFiles.add(checkDate.getDateStr());
            }
            checkDate.addDay();
        }
        return missingHydroFiles;
    }


//    /**
//     * Take a snapshot of the number of mature particles in each cell
//     * @param particles
//     * @param rp
//     * @param nSourceSites
//     * @return 
//     */
//    public static double[][] pstepMatureSnapshot(List<Particle> particles, RunProperties rp,
//            int nSourceSites) {   
//        int pstepCols = 2; 
//        // Alternatively, make a column for each source site
//        if (rp.splitPsteps == true){
//            pstepCols = nSourceSites + 1;            
//        }
//        double[][] pstepsInstMature = new double[rp.N][pstepCols];
//        for (int i = 0; i < rp.N; i++) {
//            pstepsInstMature[i][0] = i;
//        }
//        
//        for (Particle p : particles) {
//            if (p.getViable() == true) {
//                double d = 1;
//                if (rp.pstepsIncMort == true) {
//                    d = p.getDensity();
//                }
//                //System.out.println("density = "+d+" mortRate = "+p.getMortRate());
//                int elemPart = p.getElem();
//                if (rp.splitPsteps == false) {
//                    pstepsInstMature[elemPart][1] += d;//*1.0/rp.stepsPerStep;
//                } else {
//                    pstepsInstMature[elemPart][p.getStartID() + 1] += d;//*1.0/rp.stepsPerStep;
//                }
//            }
//        }
//        return pstepsInstMature;
//    } 

//    /**
//     * Take a snapshot of the number of immature particles in each cell
//     * @param particles
//     * @param rp
//     * @param nSourceSites
//     * @return 
//     */
//    public static double[][] pstepImmatureSnapshot(List<Particle> particles, RunProperties rp,
//            int nSourceSites) {   
//        int pstepCols = 2; 
//        // Alternatively, make a column for each source site
//        if (rp.splitPsteps == true){
//            pstepCols = nSourceSites + 1;
//        }
//        double[][] pstepsInstImmature = new double[rp.N][pstepCols];
//        for (int i = 0; i < rp.N; i++) {
//            pstepsInstImmature[i][0] = i;
//        }
//        
//        for (Particle p : particles) {
//            if (p.getViable() == false && p.getFree() == true) {
//                double d = 1;
//                if (rp.pstepsIncMort == true) {
//                    d = p.getDensity();
//                }
//                //System.out.println("density = "+d+" mortRate = "+p.getMortRate());
//                int elemPart = p.getElem();
//                //System.out.println("Printing to pstepsImmature");
//                if (rp.splitPsteps == false) {
//                    pstepsInstImmature[elemPart][1] += d;//*1.0/rp.stepsPerStep;
//                } else {
//                    pstepsInstImmature[elemPart][p.getStartID() + 1] += d;//*1.0/rp.stepsPerStep;
//                }
//            }
//        }
//        return pstepsInstImmature;
//    }


}
