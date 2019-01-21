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

import ucar.ma2.ArrayDouble;
import ucar.ma2.ArrayFloat;
import ucar.ma2.ArrayInt;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;
        
// import netcdf reading libraries for alternative data import
// https://publicwiki.deltares.nl/display/FEWSDOC/Reading+and+writing+netcdf+gridded+data+with+Java
//import ucar.nc2.*;

/**
 *
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
    public static int countLines(File aFile) throws Exception 
    {
        LineNumberReader reader = null;
        try {
            reader = new LineNumberReader(new FileReader(aFile));
            while ((reader.readLine()) != null);
            return reader.getLineNumber();
        } catch (Exception ex) {
            return -1;
        } finally { 
            if(reader != null) 
                reader.close();
        }
    }
    
    public static int countWords(String s){

        String trim = s.trim();
        if (trim.isEmpty()){
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
    
    public static List<HabitatSite> createHabitatSites(String filename, String sitedir, int scaleCol, boolean makeLocalCopy, List<Mesh> meshes)
    {
        List<HabitatSite> habitat = new ArrayList<>();
        System.out.println("Habitat defined in file: "+filename);
        File file = new File(sitedir+filename);
        int nLines = 0;
        try
        {
            nLines = countLines(file);
            System.out.println("FILE "+file+" NLINES "+nLines);
        }
        catch (Exception e)
        {
            System.err.println("Cannot open "+sitedir+filename);
        }
        
        // try to create habitat from file
        try
        {
            BufferedReader in = new BufferedReader(new FileReader(file));	//reading files in specified directory
 
            String line;
            //boolean printWarning=true;
            int count = 0;
            
            while ((line = in.readLine()) != null)	//file reading
            {
                //int numEntries = countWords(line);
                String[] values = line.split("\t");
                HabitatSite site = new HabitatSite(values[0],
                        (float)Double.parseDouble(values[1]),
                        (float)Double.parseDouble(values[2]),
                        (float)Double.parseDouble(values[3]),
                        (float)Double.parseDouble(values[scaleCol]),
                        meshes);
                habitat.add(site);
                count++;
            }
            in.close();
            System.out.println("Created "+count+" habitat sites");
        }
        catch (Exception e)
        {
            System.err.println("Cannot create habitat sites from "+sitedir+filename);
        }
        
        // Make a copy if required
        if (makeLocalCopy == true)
        {
            String outName = "startlocs.dat";
            try
            {
                Files.copy(file.toPath(), Paths.get(outName));
            }
            catch (Exception e)
            {
                System.err.println("Cannot copy "+sitedir+filename+" to "+outName);
            }
        }
        
        return habitat;
    }
    
    
    /**
     * Add some extra locations at which settlement is possible, or limit settlement to
     * a smaller selection
     * @param habitat
     * @param sitedir
     * @param startlocs
     * @param limit  // the maximum index of "startlocs" to allow as an endlocation
     * @return 
     */
    public static double[][] setupEndLocs(String habitat, String sitedir, double[][] startlocs, int limit)
    {
        double[][] endlocs = new double[10][3];
        
        if (habitat.equalsIgnoreCase("fishfarm3_new"))
        {
            double[][] extraEndLocs = IOUtils.readFileDoubleArray(sitedir+"160119_fishfarms_jun15.dat",205,5," ",true);
            endlocs = new double[startlocs.length+extraEndLocs.length][];
            for(int i = 0; i < startlocs.length; i++)
            {
                endlocs[i] = startlocs[i].clone();
            }
            for (int i = startlocs.length; i < startlocs.length+endlocs.length; i++)
            {
                endlocs[i] = extraEndLocs[i-startlocs.length].clone();
            }
        }
        else if (limit==0)
        {
            endlocs = new double[startlocs.length][];
            for(int i = 0; i < startlocs.length; i++)
            {
                endlocs[i] = startlocs[i].clone();
            }
        }
        else
        {
            endlocs = new double[limit][];
            for(int i = 0; i < limit; i++)
            {
                endlocs[i] = startlocs[i].clone();
            }
        }
        System.out.println("Endlocs NLINES "+endlocs.length);
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
    
    public static int[] readFileInt1D(String fullFileName) throws Exception
    {
        //double open_BC_locs[][] = new double[10][3];
        int nLines = countLines(new File(fullFileName));
        int[][] vals = IOUtils.readFileIntArray(fullFileName,nLines,1," ",true);
        int[] valsOut = new int[nLines];
        for (int i = 0; i < nLines; i++)
        {
            valsOut[i] = vals[i][0];
        }
        return valsOut;
    }
    
    
    public static float[] readNetcdfFloat1D(String filename, String variable)
    {
        System.out.println("Reading variable: "+variable);
        float[] floatOut = new float[10];
        try (NetcdfFile dataFile = NetcdfFile.open(filename, null)) {
            Variable dataVar = dataFile.findVariable(variable);
            if (dataVar == null) {
                System.out.println("Cant find Variable: "+variable);
            }   
            
            int[] shape = dataVar.getShape();
            //int origin = 0;
            
            try
            {
                ArrayFloat.D1 dataArray = (ArrayFloat.D1) dataVar.read();
                // Put the values into a native array
                floatOut = new float[shape[0]];
                for (int d1 = 0; d1 < shape[0]; d1++) {
                    floatOut[d1] = dataArray.get(d1);
                }
            }
            catch (ClassCastException e)
            {
                ArrayDouble.D1 dataArray = (ArrayDouble.D1) dataVar.read();
                // Put the values into a native array
                floatOut = new float[shape[0]];
                for (int d1 = 0; d1 < shape[0]; d1++) {
                    floatOut[d1] = (float)dataArray.get(d1);
                }
            }
            
            
        } catch (IOException ioe) {
	    ioe.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        } 
        return floatOut;
    }
    public static int[] readNetcdfInteger1D(String filename, String variable)
    {
        System.out.println("Reading variable: "+variable);
        int[] intOut = new int[10];
        try (NetcdfFile dataFile = NetcdfFile.open(filename, null)) {
            Variable dataVar = dataFile.findVariable(variable);
            if (dataVar == null) {
                System.out.println("Can't find Variable: "+variable);
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
    public static float[][] readNetcdfFloat2D(String filename, String variable, int[] origin, int[] shape)
    {
        System.out.println("Reading variable: "+variable);
        float[][] floatOut = new float[2][10];
        try (NetcdfFile dataFile = NetcdfFile.open(filename, null)) {
            Variable dataVar = dataFile.findVariable(variable);
            if (dataVar == null) {
                System.err.println("Can't find Variable: "+variable);
            }     
//            int[] origin = new int[2];
//            int[] shape = dataVar.getShape();
            
            // Check that origin is within the array
            if (origin==null 
                    || origin[0]>dataVar.getShape()[0] 
                    || origin[1]>dataVar.getShape()[1])
            {
                System.out.println("Origin not supplied, or outside bounds of variable "+variable);
                origin = new int[3];
            }
            
            //int[] shape = dataVar.getShape();
            // Check that origin+shape is within the array
            if (shape==null 
                    || origin[0]+shape[0]>dataVar.getShape()[0] 
                    || origin[1]+shape[1]>dataVar.getShape()[1])
            {
                System.out.println("Shape not supplied, or outside bounds of variable "+variable);
                shape = dataVar.getShape();
            }
            
            System.out.println(variable+" ("+shape[0]+","+shape[1]+")");
            
            try
            {
                ArrayFloat.D2 dataArray = (ArrayFloat.D2) dataVar.read(origin, shape);
                // Put the values into a native array
                floatOut = new float[shape[0]][shape[1]];
                for (int d1 = 0; d1 < shape[0]; d1++) {
                    for (int d2 = 0; d2 < shape[1]; d2++) {
                        floatOut[d1][d2] = dataArray.get(d1, d2);
                    }
                }
            }
            catch (ClassCastException e)
            {
                ArrayDouble.D2 dataArray = (ArrayDouble.D2) dataVar.read(origin, shape);
                // Put the values into a native array
                floatOut = new float[shape[0]][shape[1]];
                for (int d1 = 0; d1 < shape[0]; d1++) {
                    for (int d2 = 0; d2 < shape[1]; d2++) {
                        floatOut[d1][d2] = (float)dataArray.get(d1, d2);
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
    
    public static int[][] readNetcdfInteger2D(String filename, String variable)
    {
        System.out.println("Reading variable: "+variable);
        int[][] intOut = new int[2][10];
        try (NetcdfFile dataFile = NetcdfFile.open(filename, null)) {
            Variable dataVar = dataFile.findVariable(variable);
            if (dataVar == null) {
                System.out.println("Can't find Variable: "+variable);
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
        if(origM*origN != m*n){
            throw new IllegalArgumentException("New matrix must be of same area as matix A");
        }
        float[][] B = new float[m][n];
        float[] A1D = new float[A.length * A[0].length];

        int index = 0;
        for(int i = 0;i<A.length;i++){
            for(int j = 0;j<A[0].length;j++){
                A1D[index++] = A[i][j];
            }
        }

        index = 0;
        for(int i = 0;i<n;i++){
            for(int j = 0;j<m;j++){
                B[j][i] = A1D[index++];
            }

        }
        return B;
    }
    
    /**
     * 
     * @param filename
     * @param variable
     * @param origin
     * @param shape
     * @return
     * @throws IOException
     * @throws InvalidRangeException 
     */
    public static float[][][] readNetcdfFloat3D(String filename, String variable, int[] origin, int[] shape)
    {
        System.out.println("Reading variable: "+variable);
        float[][][] floatOut = new float[2][2][10];
        try (NetcdfFile dataFile = NetcdfFile.open(filename, null)) {
            Variable dataVar = dataFile.findVariable(variable);
            if (dataVar == null) {
                System.err.println("Can't find Variable: "+variable);
            }   
            
            // set origin externally, always?
            // set shape externally. based on known size of mesh required. Handle case when it doesn't work here, rather than needing to set here?
            
            //int[] origin = new int[3];
            // Check that origin is within the array
            if (origin==null 
                    || origin[0]>dataVar.getShape()[0] 
                    || origin[1]>dataVar.getShape()[1] 
                    || origin[2]>dataVar.getShape()[2])
            {
                System.out.println("Origin not supplied, or outside bounds of variable "+variable);
                origin = new int[3];
            }
            
            //int[] shape = dataVar.getShape();
            // Check that origin+shape is within the array
            if (shape==null 
                    || origin[0]+shape[0]>dataVar.getShape()[0] 
                    || origin[1]+shape[1]>dataVar.getShape()[1] 
                    || origin[2]+shape[2]>dataVar.getShape()[2])
            {
                System.out.println("Shape not supplied, or outside bounds of variable "+variable);
                shape = dataVar.getShape();
            }
            
            System.out.println(variable+" ("+shape[0]+","+shape[1]+","+shape[2]+")");
            
            try
            {
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
            }
            catch (ClassCastException e)
            {
                ArrayDouble.D3 dataArray = (ArrayDouble.D3) dataVar.read(origin, shape);
                // Put the values into a native array
                floatOut = new float[shape[0]][shape[1]][shape[2]];
                for (int d1 = 0; d1 < shape[0]; d1++) {
                    for (int d2 = 0; d2 < shape[1]; d2++) {
                        for (int d3 = 0; d3 < shape[2]; d3++) {
                            floatOut[d1][d2][d3] = (float)dataArray.get(d1, d2, d3);
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
     * 
     * @param filename
     * @param variable
     * @param origin
     * @param shape
     * @return
     * @throws IOException
     * @throws InvalidRangeException 
     */
    public static float[][][][] readNetcdfFloat4D(String filename, String variable, int[] origin, int[] shape)
    {
        System.out.println("Reading variable: "+variable);
        float[][][][] floatOut = new float[2][2][10][1];
        try (NetcdfFile dataFile = NetcdfFile.open(filename, null)) {
            Variable dataVar = dataFile.findVariable(variable);
            if (dataVar == null) {
                System.err.println("Can't find Variable: "+variable);
            }   
            
            // set origin externally, always?
            // set shape externally. based on known size of mesh required. Handle case when it doesn't work here, rather than needing to set here?
            
            //int[] origin = new int[3];
            // Check that origin is within the array
            if (origin==null 
                    || origin[0]>dataVar.getShape()[0] 
                    || origin[1]>dataVar.getShape()[1] 
                    || origin[2]>dataVar.getShape()[2]
                    || origin[3]>dataVar.getShape()[3])
            {
                System.out.println("Origin not supplied, or outside bounds of variable "+variable);
                origin = new int[3];
            }
            
            //int[] shape = dataVar.getShape();
            // Check that origin+shape is within the array
            if (shape==null 
                    || origin[0]+shape[0]>dataVar.getShape()[0] 
                    || origin[1]+shape[1]>dataVar.getShape()[1] 
                    || origin[2]+shape[2]>dataVar.getShape()[2]
                    || origin[3]+shape[3]>dataVar.getShape()[3])
            {
                System.out.println("Shape not supplied, or outside bounds of variable "+variable);
                shape = dataVar.getShape();
            }
            
            System.out.println(variable+" ("+shape[0]+","+shape[1]+","+shape[2]+","+shape[3]+")");
            
            try
            {
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
            }
            catch (ClassCastException e)
            {
                ArrayDouble.D4 dataArray = (ArrayDouble.D4) dataVar.read(origin, shape);
                // Put the values into a native array
                floatOut = new float[shape[0]][shape[1]][shape[2]][shape[3]];
                for (int d1 = 0; d1 < shape[0]; d1++) {
                    for (int d2 = 0; d2 < shape[1]; d2++) {
                        for (int d3 = 0; d3 < shape[2]; d3++) {
                            for (int d4 = 0; d4 < shape[2]; d4++) {
                                floatOut[d1][d2][d3][d3] = (float)dataArray.get(d1, d2, d3, d3);
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
    
    
    public static double[][] readFileDoubleArray(String filename, int rows, int cols, String sep, boolean note)
    {
        double[][] myDouble = new double[rows][cols];
        int x=0, y=0;
        boolean failed = false;
        double sum = 0;
        try
        {
            BufferedReader in = new BufferedReader(new FileReader(filename));	//reading files in specified directory
 
            String line;
            boolean printWarning=true;
            
            outer: while ((line = in.readLine()) != null)	//file reading
            {
                
                int numEntries = countWords(line);
                if (numEntries < cols && printWarning==true)
                {
                    System.out.println("WARNING: Number of entries on line = "+numEntries);
                    printWarning=false;
                }
                y=0;
                if (x >= rows)
                {
                    System.out.println(filename+" has more rows than expected.");
                    failed = true;
                    break;
                }
                String[] values = line.split(" ");
                for (String str : values)
                {
                    if (y >= cols)
                    {
                        System.out.println(filename+" has more columns than expected: "+y);
                        failed = true;
                        break outer;
                    }
                    double str_double = Double.parseDouble(str);
                    myDouble[x][y]=str_double;
                    sum+=myDouble[x][y];
                    //System.out.print(myDouble[x][y] + " ");
                    y++;

                }

                //System.out.println("");
                x++;
            }
            in.close();
            
        } catch( IOException ioException ) {
            throw new RuntimeException(ioException);
            //System.err.println("******************* Cannot read from file "+filename+" ******************************");
            //failed = true;
        }
        if (note == true && failed == false)
        {
            System.out.printf("Created %dx%d array from file: %s\n",myDouble.length,myDouble[0].length,filename);
            System.out.println("Array sum at read time = "+sum);
        }
        else if (failed == true)
        {
            System.out.println("FAILED to read file "+filename);
            System.exit(1);
        }
        return myDouble;
    }
    
    public static int[][] readFileIntArray(String filename, int rows, int cols, String sep, boolean note)
    {
        int[][] myInt = new int[rows][cols];
        int x=0, y=0;
        boolean failed = false;
        try
        {
            BufferedReader in = new BufferedReader(new FileReader(filename));	//reading files in specified directory
            //System.out.println("in readFileIntArray");
            String line;
            while ((line = in.readLine()) != null)	//file reading
            {
                y=0;
                if (x >= rows)
                {
                    System.out.println("WARNING: "+filename+" has more rows than expected. EXITING!");
                    failed = true;
                    break;
                }
                String[] values = line.split(sep);
                for (String str : values)
                {
                    if (y >= cols)
                    {
                        System.out.println("WARNING: "+filename+" has more columns than expected. Some information not read!");
                        //failed = true;
                        //break;
                    }
                    int str_int = Integer.parseInt(str);
                    myInt[x][y]=str_int;
                    //System.out.print(myDouble[x][y] + " ");
                    y++;

                }
                //System.out.println("");
                x++;
            }
            in.close();
            
        
        } catch( Exception e ) {
            System.err.println("******************* Cannot read from file "+filename+" ******************************");
            failed = true;
        }
        if (note == true && failed == false)
        {
            System.out.printf("Created %dx%d array from file: %s\n",myInt.length,myInt[0].length,filename);
        }
        else if (failed == true)
        {
            System.out.println("FAILED to read file "+filename);
            System.exit(1);
        }
        return myInt;
    }
    
    /**
     * Print a predetermined string of characters to a file. Also ensures a new 
     * file is started for a given day for particle locations
     * @param headerString
     * @param filename 
     */
    public static void printFileHeader(String headerString, String filename)
    {
        try
        {
            FileWriter fstream = new FileWriter(filename);
            PrintWriter out = new PrintWriter(fstream);
            out.printf(headerString);
            out.printf("\n");
            out.close();
        } 
        catch (Exception e)
        {//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }
    
    /**
     * Print an array of doubles directly to a file, 
     * @param variable
     * @param filename 
     * @param asInt
     */
    public static void writeFloatArrayToFile(float[][] variable, String filename, boolean asInt)
    {
        try
        {
            // Create file 
            FileWriter fstream = new FileWriter(filename);
            PrintWriter out = new PrintWriter(fstream);
            for (int i = 0; i < variable.length; i++)
            {
                for (int j = 0; j < variable[0].length; j++)
                {   
                    if (asInt==false)
                    {
                        out.printf("%.4e ",variable[i][j]);
                    }
                    else
                    {
                        out.printf("%d ",(int)variable[i][j]);
                    }
                }
                out.printf("\n");
            }
            //Close the output stream
            out.close();
        }catch (Exception e){//Catch exception if any
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
    public static void writeIntegerArrayToFile(int[][] variable, String filename)
    {
        try
        {
            // Create file 
            FileWriter fstream = new FileWriter(filename);
            PrintWriter out = new PrintWriter(fstream);
            for (int i = 0; i < variable.length; i++)
            {
                for (int j = 0; j < variable[0].length; j++)
                {
                    out.printf("%d ",variable[i][j]);
                }
                out.printf("\n");
            }
            //Close the output stream
            out.close();
        }catch (Exception e){//Catch exception if any
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
    public static void particleLocsToFile(Particle[] particles, int npartsSaved, int tt, String filename)
    {
        try
        {
            // Create file 
            FileWriter fstream = new FileWriter(filename,true);
            PrintWriter out = new PrintWriter(fstream);
            for (int i = 0; i < npartsSaved; i++)
            {
                out.printf("%d %d %.1f %.1f %d %d\n",tt,particles[i].getID(),
                        particles[i].getLocation()[0],particles[i].getLocation()[1],
                        particles[i].getElem(),particles[i].getStatus());
            }
            //Close the output stream
            out.close();
        }catch (Exception e){//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }
    
    public static void particleLocsToFile(List<Particle> particles, int npartsSaved, int tt, String filename)
    {
        try
        {
            // Create file 
            FileWriter fstream = new FileWriter(filename,true);
            PrintWriter out = new PrintWriter(fstream);
            for (int i = 0; i < npartsSaved; i++)
            {
                out.printf("%d %d %.1f %.1f %d %d\n",tt,particles.get(i).getID(),
                        particles.get(i).getLocation()[0],particles.get(i).getLocation()[1],
                        particles.get(i).getElem(),particles.get(i).getStatus());
            }
            //Close the output stream
            out.close();
        }catch (Exception e){//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }
    
    /** Print ALL particle locations at a particular time, with corresponding start locations.
     * 
     * @param particles
     * @param filename 
     */
    public static void particleLocsToFile1(Particle[] particles, String filename)
    {
        try
        {
            // Create file 
            FileWriter fstream = new FileWriter(filename,true);
            PrintWriter out = new PrintWriter(fstream);
            for (int i = 0; i < particles.length; i++)
            {
                out.printf("%d %f %f %f %f\n",i,particles[i].getStartLocation()[0],particles[i].getStartLocation()[1],
                        particles[i].getLocation()[0],particles[i].getLocation()[1]);
            }
            //Close the output stream
            out.close();
        }catch (Exception e){//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }
    public static void particleLocsToFile1(List<Particle> particles, String filename, boolean append)
    {
        try
        {
            // Create file 
            FileWriter fstream = new FileWriter(filename,append);
            PrintWriter out = new PrintWriter(fstream);
            for (int i = 0; i < particles.size(); i++)
            {
                out.printf("%d %f %f %f %f\n",i,particles.get(i).getStartLocation()[0],
                        particles.get(i).getStartLocation()[1],particles.get(i).getLocation()[0],particles.get(i).getLocation()[1]);
            }
            //Close the output stream
            out.close();
        }catch (Exception e){//Catch exception if any
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
    public static void particleLocsToFile_full(List<Particle> particles, int currentHour, String filename, boolean append)
    {
        try
        {
            // Create file 
            FileWriter fstream = new FileWriter(filename,append);
            PrintWriter out = new PrintWriter(fstream);
            for (Particle p : particles)
            {
                out.printf("%d %d %s %.1f %s %.4f %.4f %d %d %.2f\n",
                        currentHour,
                        p.getID(),
                        p.getStartDate().getDateStr(),
                        p.getAge(),
                        p.getStartID(),
                        p.getLocation()[0],
                        p.getLocation()[1],
                        p.getElem(),
                        p.getStatus(),
                        p.getDensity()
                );
            }
            
            //Close the output stream
            out.close();
        }catch (Exception e){//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }
    
    
    public static void particleLocsToNetcdfFile()
    {
        // https://www.unidata.ucar.edu/software/netcdf/examples/programs/Simple_xy_wr.java
        // What about appending data? Need to check how that would work.
    }
    
    
    
    public static void arrivalToFile(Particle p, ISO_datestr currentDate, double currentTime, String filename, boolean append)
    {
        try
        {
            // Create file 
            FileWriter fstream = new FileWriter(filename,append);
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
        }catch (Exception e){//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }
    
    
    /**
     * Print (append) single particle location to file
     * @param time
     * @param particles
     * @param i
     * @param filename 
     */
    public static void particleLocsToFile2(double time, Particle[] particles, int i, String filename)
    {
        try
        {
            // Create file 
            FileWriter fstream = new FileWriter(filename,true);
            PrintWriter out = new PrintWriter(fstream);
            out.printf("%f %f %f %d\n",time,particles[i].getLocation()[0],particles[i].getLocation()[1],particles[i].getElem());
            //Close the output stream
            out.close();
        }catch (Exception e){//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }
    
    
    //    /**
//     * Make additions to the element presence counts (PSTEPS)
//     *
//     * @param particles
//     * @param rp
//     * @param pstepsMature
//     * @param pstepsImmature
//     * @param subStepDt
//     */
//    public static void pstepUpdater(List<Particle> particles, RunProperties rp,
//            double[][] pstepsMature, double[][] pstepsImmature, double subStepDt) {
//        for (Particle p : particles) {
//            double d = 1;
//            if (rp.pstepsIncMort == true) {
//                d = p.getDensity();
//            }
//            //System.out.println("density = "+d+" mortRate = "+p.getMortRate());
//            int elemPart = p.getElem();
//            // psteps arrays are updated by lots of threads
//            if (p.getViable() == true) {
//                if (rp.splitPsteps == false) {
//                    pstepsMature[elemPart][1] += d * (subStepDt / 3600);//*1.0/rp.stepsPerStep;
//                } else {
//                    pstepsMature[elemPart][p.getStartID() + 1] += d * (subStepDt / 3600);//*1.0/rp.stepsPerStep;
//                }
//            } else if (p.getFree() == true) {
//                //System.out.println("Printing to pstepsImmature");
//                if (rp.splitPsteps == false) {
//                    pstepsImmature[elemPart][1] += d * (subStepDt / 3600);//*1.0/rp.stepsPerStep;
//                } else {
//                    pstepsImmature[elemPart][p.getStartID() + 1] += d * (subStepDt / 3600);//*1.0/rp.stepsPerStep;
//                }
//            }
//        }
//    }
    
    /**
     * Take a snapshot of the number of mature particles in each cell
     * @param particles
     * @param rp
     * @param nSourceSites
     * @return 
     */
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
