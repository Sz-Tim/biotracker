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

import java.util.*;
        
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
    
    public static double[][] setupStartLocs(String filename, String sitedir)
    {
        double startlocs[][] = new double[10][3];
        System.out.println("Startlocations defined in file: "+filename);
        File file = new File(sitedir+filename);
        int nLines = 0;
        try
        {
            nLines = countLines(file);
            System.out.println("FILE "+file+" NLINES "+nLines);
        }
        catch (Exception e)
        {
            System.out.println("Cannot open "+sitedir+filename);
        }
        startlocs = IOUtils.readFileDoubleArray(sitedir+filename,nLines,5," ",true);
        return startlocs;
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
    
    public static double[][] setupOpenBCLocs(String location, String datadir2)
    {
        double open_BC_locs[][] = new double[10][3];
        if (location.equalsIgnoreCase("minch") || location.equalsIgnoreCase("minch_continuous") || location.equalsIgnoreCase("minch_jelly"))
        {
            //String sitedir=basedir+"minch_sites/";
            //open_BC_locs = IOUtils.readFileDoubleArray(sitedir+"open_boundary_locs.dat",120,3," ",true);
            open_BC_locs = IOUtils.readFileDoubleArray(datadir2+"open_boundary_locs_os.dat",137,3," ",true);
        } 
        else
        {
            open_BC_locs = IOUtils.readFileDoubleArray("C:\\Users\\sa01ta\\Documents\\lorn\\120903_renewableimpact\\131107_revision\\open_boundary_locs.dat",93,3," ",true);
        }
        return open_BC_locs;
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
    
 
    public static void writeDoubleArrayToFile(double[][] variable, String filename)
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
                    out.printf("%.4e ",variable[i][j]);
                }
                out.printf("\n");
            }
            //Close the output stream
            out.close();
        }catch (Exception e){//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }
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
                out.printf("%d %d %f %f %d\n",tt,particles[i].getID(),particles[i].getLocation()[0],particles[i].getLocation()[1],particles[i].getElem());
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
                out.printf("%d %d %f %f %d\n",tt,particles.get(i).getID(),
                        particles.get(i).getLocation()[0],particles.get(i).getLocation()[1],particles.get(i).getElem());
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
                out.printf("%d %f %f %f %f\n",i,particles[i].getStartLocation()[0],particles[i].getStartLocation()[1],particles[i].getLocation()[0],particles[i].getLocation()[1]);
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
}
