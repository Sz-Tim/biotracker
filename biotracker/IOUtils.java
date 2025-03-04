/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotracker;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.io.LineNumberReader;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;

import java.util.*;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.filefilter.DirectoryFileFilter;
import org.apache.commons.io.filefilter.RegexFileFilter;
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
@SuppressWarnings("ConstantConditions")
public class IOUtils {

    /**
     * Count the lines in a file, in order to read files of unknown length
     */
    @SuppressWarnings("StatementWithEmptyBody")
    public static int countLines(File aFile) {
        try (LineNumberReader reader = new LineNumberReader(new FileReader(aFile))) {
            while ((reader.readLine()) != null) ;
            return reader.getLineNumber();
        } catch (Exception ex) {
            return -1;
        }
    }

    public static int countWords(String s) {
        String trim = s.trim();
        if (trim.isEmpty()) {
            return 0;
        }
        return trim.split("\\w+").length;
    }

    public static List<HabitatSite> createHabitatSites(String filename, String sitedir, int scaleCol, boolean makeLocalCopy, List<Mesh> meshes, RunProperties rp) {
        List<HabitatSite> habitat = new ArrayList<>();
        System.out.println("Habitat defined in file: " + filename);
        File file;

        if (sitedir != null) {
            file = new File(sitedir + filename);
        } else {
            file = new File(filename);
        }

        int nLines;
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
            boolean skipHeader = true;
            int count = 0, countNotCreated = 0;

            while ((line = in.readLine()) != null) {
                if(skipHeader) {
                    skipHeader = false;
                    continue;
                }
                String[] values = line.split(",");
                String ID = values[0];
                float x = (float) Double.parseDouble(values[1]);
                float y = (float) Double.parseDouble(values[2]);
                float depth = 0;
                if (values.length > 3) {
                    depth = (float) Double.parseDouble(values[3]);
                }
                float scale = 1;
                if (values.length > 4 && scaleCol < values.length) {
                    scale = (float) Double.parseDouble(values[scaleCol]);
                }
                HabitatSite site = new HabitatSite(ID, x, y, depth, scale, meshes, rp);
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
     * @param limit     // the maximum index of "startlocs" to allow as an endlocation
     */
    public static double[][] setupEndLocs(String habitat, String sitedir, double[][] startlocs, int limit) {
        double[][] endlocs;

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

    public static List<Particle> readRestartParticles(RunProperties rp) {
        List<Particle> parts = new ArrayList<>();
        System.out.println("Particles to restart defined in file: " + rp.restartParticles);
        File file = new File(rp.restartParticles);
        int nLines;
        try {
            nLines = countLines(file);
            System.out.println("FILE " + file + " NLINES " + nLines);
        } catch (Exception e) {
            System.err.println("Cannot open " + rp.restartParticles);
        }

        int[] todayInt = ISO_datestr.dateIntParse(rp.start_ymd);
        ISO_datestr today = new ISO_datestr(todayInt[0], todayInt[1], todayInt[2]);
        int todayDateNum = today.getDateNum();

        try {
            BufferedReader in = new BufferedReader(new FileReader(file));    //reading files in specified directory

            String line;
            int count = 0;

            while ((line = in.readLine()) != null)    //file reading
            {
                // Ignore the header line
                if (count > 0) {
                    Particle p = new Particle(line, rp.species);
                    boolean add = true;
                    if (rp.restartParticlesCutoffDays > 0) {
                        ISO_datestr pStart = p.getStartDate();
                        int pStartDateNum = pStart.getDateNum();
                        if (todayDateNum - pStartDateNum > rp.restartParticlesCutoffDays) {
                            add = false;
                        }
                    }
                    if (add) {
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


    public static double[][] readDailyDensities(String filename, List<String> siteStartNames, RunProperties rp) {
        double[][] dailyDensities = new double[siteStartNames.size()][rp.numberOfDays];
        int nLines;
        try {
            nLines = countLines(new File(filename));
            System.out.println("FILE " + filename + " NLINES " + nLines);
        } catch (Exception e) {
            System.err.println("Cannot open " + filename);
        }

        // try to create habitat from file
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));    //reading files in specified directory
            String line;
            boolean skipHeader = true;

            while ((line = in.readLine()) != null) {
                if(skipHeader) {
                    skipHeader = false;
                    continue;
                }
                String[] values = line.split(",");
                String ID = values[0];
                int site = siteStartNames.indexOf(ID);
                for (int i=0; i < rp.numberOfDays; i++) {
                    dailyDensities[site][i] = Float.parseFloat(values[i+1]);
                }
            }
            in.close();
        } catch (Exception e) {
            System.err.println("Cannot load daily site densities from " + filename);
        }
        return(dailyDensities);
    }


    public static int[] readFileInt1D(String fullFileName) {
        int nLines = countLines(new File(fullFileName));
        int[][] vals = IOUtils.readFileIntArray(fullFileName, nLines, 1, ",", true);
        int[] valsOut = new int[nLines];
        for (int i = 0; i < nLines; i++) {
            valsOut[i] = vals[i][0];
        }
        return valsOut;
    }


    public static double[] readFileDouble1D(String fullFileName) {
        int nLines = countLines(new File(fullFileName));
        double[][] vals = IOUtils.readFileDoubleArray(fullFileName, nLines, 1, ",", false);
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


        } catch (Exception ioe) {
            ioe.printStackTrace();
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
        } catch (Exception ioe) {
            ioe.printStackTrace();
        }
        return intOut;
    }

    /**
     * Method to read NetCDF variable
     */
    public static float[][] readNetcdfFloat2D(String filename, String variable, int[] origin, int[] shape) {
        float[][] floatOut = null;

        try (NetcdfFile dataFile = NetcdfFile.open(filename, null)) {
            Variable dataVar = dataFile.findVariable(variable);
            if (dataVar == null) {
                System.err.println("  Can't find Variable: " + variable);
                return (floatOut);
            }

            // Check that origin is within the array
            if (origin == null
                    || origin[0] > dataVar.getShape()[0]
                    || origin[1] > dataVar.getShape()[1]) {
                origin = new int[3];
            }

            // Check that origin+shape is within the array
            if (shape == null
                    || origin[0] + shape[0] > dataVar.getShape()[0]
                    || origin[1] + shape[1] > dataVar.getShape()[1]) {
                shape = dataVar.getShape();
            }

            System.out.println(variable + " (" + shape[0] + "," + shape[1] + ")");

            try {
                ArrayFloat.D2 dataArray = (ArrayFloat.D2) dataVar.read(origin, shape);
                // Put the values into a native array
                floatOut = new float[shape[0]][shape[1]];
                floatOut = (float[][]) dataArray.copyToNDJavaArray();
            } catch (ClassCastException e) {
                ArrayDouble.D2 dataArray = (ArrayDouble.D2) dataVar.read(origin, shape);
                // Put the values into a native array
                floatOut = new float[shape[0]][shape[1]];
                floatOut = (float[][]) dataArray.copyToNDJavaArray();
            }

        } catch (Exception ioe) {
            ioe.printStackTrace();
        }
        return floatOut;
    }

    public static int[][] readNetcdfInteger2D(String filename, String variable) {
        System.out.println("Reading variable: " + variable);
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
        } catch (Exception ioe) {
            ioe.printStackTrace();
        }
        return intOut;
    }

    /**
     * Reshape a two-dimensional array of floating point numbers
     */
    public static float[][] reshapeFloat(float[][] A, int m, int n) {
        int origM = A.length;
        int origN = A[0].length;
        if (origM * origN != m * n) {
            throw new IllegalArgumentException("New matrix must be of same area as matrix A");
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


    public static float[][][] readNetcdfFloat3D(String filename, String variable, int[] origin, int[] shape) {
        float[][][] floatOut = null;

        try (NetcdfFile dataFile = NetcdfFile.open(filename, null)) {
            Variable dataVar = dataFile.findVariable(variable);
            if (dataVar == null) {
                System.err.println("  Can't find Variable: " + variable);
                return (floatOut);
            }

            // Check that origin is within the array
            if (origin == null
                    || origin[0] > dataVar.getShape()[0]
                    || origin[1] > dataVar.getShape()[1]
                    || origin[2] > dataVar.getShape()[2]) {
                origin = new int[3];
            }

            // Check that origin+shape is within the array
            if (shape == null
                    || origin[0] + shape[0] > dataVar.getShape()[0]
                    || origin[1] + shape[1] > dataVar.getShape()[1]
                    || origin[2] + shape[2] > dataVar.getShape()[2]) {
                shape = dataVar.getShape();
            }

            System.out.println(variable + " (" + shape[0] + "," + shape[1] + "," + shape[2] + ")");

            try {
                ArrayFloat.D3 dataArray = (ArrayFloat.D3) dataVar.read(origin, shape);
                floatOut = new float[shape[0]][shape[1]][shape[2]];
                floatOut = (float[][][]) dataArray.copyToNDJavaArray();
            } catch (ClassCastException e) {
                ArrayDouble.D3 dataArray = (ArrayDouble.D3) dataVar.read(origin, shape);
                // Put the values into a native array
                floatOut = new float[shape[0]][shape[1]][shape[2]];
                floatOut = (float[][][]) dataArray.copyToNDJavaArray();
            }

        } catch (Exception ioe) {
            ioe.printStackTrace();
        }
        return floatOut;
    }


    public static float[][][][] readNetcdfFloat4D(String filename, String variable, int[] origin, int[] shape) {
        System.out.println("Reading variable: " + variable);
        float[][][][] floatOut = null;

        try (NetcdfFile dataFile = NetcdfFile.open(filename, null)) {
            Variable dataVar = dataFile.findVariable(variable);
            if (dataVar == null) {
                System.err.println("  Can't find Variable: " + variable);
                return (floatOut);
            }
            // Check that origin is within the array
            if (origin == null
                    || origin[0] > dataVar.getShape()[0]
                    || origin[1] > dataVar.getShape()[1]
                    || origin[2] > dataVar.getShape()[2]
                    || origin[3] > dataVar.getShape()[3]) {
                System.out.println("Origin not supplied, or outside bounds of variable " + variable);
                origin = new int[3];
            }

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

        } catch (Exception ioe) {
            ioe.printStackTrace();
        }
        return floatOut;
    }


    public static double[][] readFileDoubleArray(String filename, int rows, int cols, String sep, boolean note) {
        double[][] myDouble = new double[rows][cols];
        int x = 0, y;
        boolean failed = false;
        double sum = 0;
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));    //reading files in specified directory

            String line;
            boolean printWarning = true;

            outer:
            while ((line = in.readLine()) != null) {
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
                    y++;

                }
                x++;
            }
            in.close();

        } catch (IOException ioException) {
            throw new RuntimeException(ioException);
        }
        if (note && !failed) {
            System.out.printf("Created %dx%d array from file: %s\n", myDouble.length, myDouble[0].length, filename);
            System.out.println("Array sum at read time = " + sum);
        } else if (failed) {
            System.err.println("FAILED to read file " + filename);
            System.exit(1);
        }
        return myDouble;
    }

    public static int[][] readFileIntArray(String filename, int rows, int cols, String sep, boolean note) {
        int[][] myInt = new int[rows][cols];
        int x = 0, y;
        boolean failed = false;
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));    //reading files in specified directory
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
                    }
                    int str_int = Integer.parseInt(str);
                    myInt[x][y] = str_int;
                    y++;
                }
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
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));    //reading files in specified directory

            String line;
            boolean printWarning = true;

            while ((line = in.readLine()) != null) {
                int numEntries = countWords(line);
                if (numEntries < 3 && printWarning) {
                    System.out.println("WARNING: Number of entries on line = " + numEntries);
                    printWarning = false;
                }
                String[] values = line.split(",");
                ISO_datestr rowDate = new ISO_datestr(values[0]);
                if (rowDate.getDateNum() >= startDate.getDateNum() && rowDate.getDateNum() <= endDate.getDateNum()) {
                    if (note) {
                        System.out.println("Adding sunrise & sunset for " + rowDate + ": " + values[1] + " " + values[2]);
                    }
                    daylightHours[dayCount][0] = Integer.parseInt(values[1]);
                    daylightHours[dayCount][1] = Integer.parseInt(values[2]);
                    dayCount++;
                }

            }
            in.close();
        } catch (IOException ioException) {
            throw new RuntimeException(ioException);
        }
        if (note) {
            System.out.printf("Created %dx%d array from file: %s\n", daylightHours.length, daylightHours[0].length, filename);
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
    // ID date hour step startLocation startDate age density x y z layer status degreeDays sink swim temp salinity mortality dX dY dZ
    public static void writeMovements(Particle part, String currentDate, int currentHour, int step,
                                      int sink, int swim, double temperature, double temperatureSurface, double salinity, double[] displacement,
                                      String filename, boolean append) {
        try {
            // Create file
            FileWriter fstream = new FileWriter(filename, append);
            PrintWriter out = new PrintWriter(fstream);
            out.printf("%d,%s,%d,%d,%s,%.4f,%.4f,%.1f,%.1f,%.4f,%d,%d,%d,%.4f,%d,%d,%.4f,%.4f,%.4f,%.4f",
                    part.getID(),
                    currentDate,
                    currentHour,
                    step,
                    part.getStartDate().getDateStr(),
                    part.getAge(),
                    part.getDensity(),
                    part.getLocation()[0],
                    part.getLocation()[1],
                    part.getDepth(),
                    part.getDepthLayer(),
                    part.getMesh(),
                    part.getStatus(),
                    part.getDegreeDays(),
                    sink, swim,
                    temperature,
                    salinity,
                    part.getMortRate(),
                    temperatureSurface
                    );
            for (double dim: displacement) {
                out.printf(",%.7f", dim);
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
     * Print an array of doubles directly to a file
     */
    public static void writeFloatArrayToFile(float[][] variable, String filename, boolean asInt, boolean firstColInt) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(filename);
            PrintWriter out = new PrintWriter(fstream);
            for (float[] floats : variable) {
                for (int j = 0; j < variable[0].length; j++) {
                    if (asInt || (firstColInt && j == 0)) {
                        out.printf("%d", (int) floats[j]);
                    } else {
                        out.printf("%.4e", floats[j]);
                    }
                    if (j < variable[0].length - 1) {
                        out.printf(",");
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


    public static void writeIntegerArrayToFile(int[][] variable, String filename) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(filename);
            PrintWriter out = new PrintWriter(fstream);
            for (int[] ints : variable) {
                for (int j = 0; j < variable[0].length; j++) {
                    out.printf("%d,", ints[j]);
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
     * @param npartsSaved the number of particles for which to save locations (from id=0)
     */
    public static void particleLocsToFile(Particle[] particles, int npartsSaved, int tt, String filename) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(filename, true);
            PrintWriter out = new PrintWriter(fstream);
            for (int i = 0; i < npartsSaved; i++) {
                out.printf("%d,%d,%.1f,%.1f,%.1f,%d,%d,%d\n", tt, particles[i].getID(),
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
                out.printf("%d,%d,%.1f,%.1f,%.1f,%d,%d,%d\n", tt, particles.get(i).getID(),
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
     */
    public static void particleLocsToFile1(Particle[] particles, String filename) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(filename, true);
            PrintWriter out = new PrintWriter(fstream);
            for (int i = 0; i < particles.length; i++) {
                out.printf("%d,%f,%f,%f,%d,%f,%f\n", i,
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
                out.printf("%d,%f,%f,%f,%f,%f,%d\n", i,
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
     */
    public static void particleLocsToFile_full(List<Particle> particles, int currentHour, String filename, boolean append) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(filename, append);
            PrintWriter out = new PrintWriter(fstream);
            for (Particle p : particles) {
                out.printf("%d,%d,%s,%.1f,%s,%.1f,%.1f,%.1f,%d,%d,%d,%.4f,%d\n",
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


    public static void particlesToRestartFile(List<Particle> particles, int currentHour, String filename, boolean append, RunProperties rp, int partSubset) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(filename, append);
            PrintWriter out = new PrintWriter(fstream);
            for (Particle p : particles) {
                if (!rp.coordOS) {
                    if (p.getID() % partSubset == 0){
                        out.printf("%d,%d,%s,%.1f,%s,%.7f,%.7f,%d,%d,%.4f,%d,%.2f,%d,%.2f,%.2f,%.2f,%.2f,%.2f\n",
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
                                p.getxyTotal(),
                                p.getzTotal()
                        );
                    }
                } else {
                    if (p.getID() % partSubset == 0) {
                        out.printf("%d,%d,%s,%.1f,%s,%.1f,%.1f,%d,%d,%.4f,%d,%.2f,%d,%.2f,%.2f,%.2f,%.2f,%.2f\n",
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
                                p.getxyTotal(),
                                p.getzTotal()
                        );
                    }
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


    /**
     * Print (append) single particle location to file
     */
    public static void particleLocsToFile2(double time, Particle[] particles, int i, String filename) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(filename, true);
            PrintWriter out = new PrintWriter(fstream);
            out.printf("%f,%f,%f,%d\n", time, particles[i].getLocation()[0], particles[i].getLocation()[1], particles[i].getElem());
            //Close the output stream
            out.close();
        } catch (Exception e) {//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }


    /**
     * Make additions to the element presence counts (PSTEPS)
     * Note that this ONLY applies to the first mesh when there are nested meshes!
     */
    public static void pstepsUpdater(List<Particle> particles, RunProperties rp,
                                     float[][] pstepsMature, float[][] pstepsImmature, double dt) {
        for (Particle p : particles) {
            if(p.getMesh() > 0) {
                continue;
            }
            if(p.getDepth() > rp.pstepsMaxDepth) {
                continue;
            }
            double d = p.getDensity();
            int elemPart = p.getElem();
            int col = rp.splitPsteps ? p.getStartIndex() : 0;
            if (p.isInfectious()) {
                pstepsMature[elemPart][col] += (float) (d * (dt / 3600));
            } else if (rp.recordImmature && p.isFree()) {
                pstepsImmature[elemPart][col] += (float) (d * (dt / 3600));
            }
        }
    }

    public static void pstepsSparseUpdater(List<Particle> particles, RunProperties rp,
                                           List<SparseFloatArray> pstepsMature, List<SparseFloatArray> pstepsImmature, double subStepDt) {
        if(rp.pstepsMaxDepth > 1000) {
            for (Particle p : particles) {
                if(p.getMesh() > 0) {
                    continue;
                }
                pstepsAddParticle(rp, pstepsMature, pstepsImmature, subStepDt, p);
            }
        } else {
            for (Particle p : particles) {
                if(p.getMesh() > 0) {
                    continue;
                }
                if (p.getDepth() > rp.pstepsMaxDepth) {
                    continue;
                }
                pstepsAddParticle(rp, pstepsMature, pstepsImmature, subStepDt, p);
            }
        }
    }

    public static void pstepsAddParticle(RunProperties rp, List<SparseFloatArray> pstepsMature, List<SparseFloatArray> pstepsImmature, double subStepDt, Particle p) {
        double d = p.getDensity();
        int elemPart = p.getElem();
        int col = rp.splitPsteps ? p.getStartIndex() + 1 : 1;
        if (p.isInfectious()) {
            SparseFloatArray arr = pstepsMature.get(col); // Get the relevant sparse array from the list
            arr.put(col, arr.get(elemPart) + (float) (d * (subStepDt / 3600))); // place the new value in the array
            pstepsMature.set(col, arr); // put the array back in the list

        } else if (rp.recordImmature && p.isFree()) {
            SparseFloatArray arr = pstepsImmature.get(col); // Get the relevant sparse array from the list
            arr.put(col, arr.get(elemPart) + (float) (d * (subStepDt / 3600))); // place the new value in the array
            pstepsImmature.set(col, arr); // put the array back in the list
        }
    }

    /**
     * Make additions to the element vertical distribution
     * NOTE: Only works with first mesh when nested meshes are used!
     */
    public static void vertDistrUpdater(List<Particle> particles, RunProperties rp,
                                        float[][] vertDistrMature, float[][] vertDistrImmature, double dt) {
        for (Particle p : particles) {
            if(p.getMesh() > 0) {
                continue;
            }
            double d = p.getDensity();
            int elemPart = p.getElem();
            int col = (int) Math.min(Math.floor(p.getDepth()), rp.vertDistrMax);
            if (p.isInfectious()) {
                vertDistrMature[elemPart][col] += (float) (d * dt / 3600);
            } else if (rp.recordImmature && p.isFree()) {
                vertDistrImmature[elemPart][col] += (float) (d * dt / 3600);
            }
        }
    }

    public static void writeNonZeros2DArrayToCSV(float[][] data, String header, String rowFormat, String filename) throws IOException {
        try (Writer writer = new FileWriter(filename)) {
            writer.write(header + "\n"); // Header row
            for (int row = 0; row < data.length; row++) {
                for (int col = 0; col < data[row].length; col++) {
                    if (data[row][col] != 0) {
                        writer.write(String.format(rowFormat + "\n", row, col, data[row][col]));
                    }
                }
            }
        }
    }


    public static ArrayList<String> checkHydroFilesExist(RunProperties rp, ISO_datestr startDate, ISO_datestr endDate, int numberOfDays) {
        ArrayList<String> missingHydroFiles = new ArrayList<>();
        String filePathPrefix = rp.mesh1Domain;

        File file1 = new File(rp.datadir + System.getProperty("file.separator"));
        List<File> allMatchingFiles = (List<File>) FileUtils.listFiles(
                file1,
                new RegexFileFilter(rp.mesh1Domain + ".*nc"),
                DirectoryFileFilter.DIRECTORY); // Optional directory filter

        ISO_datestr checkDate = new ISO_datestr(startDate.getDateStr());
        for (int i = 0; i < numberOfDays; i++) {
            String filename = String.format("%s_%04d%02d%02d", filePathPrefix, checkDate.getYear(), checkDate.getMonth(), checkDate.getDay());
            if (!allMatchingFiles.stream().anyMatch(f -> f.getName().contains(filename))) {
                missingHydroFiles.add(checkDate.getDateStr());
            }
            checkDate.addDay();
        }
        return missingHydroFiles;
    }


}
