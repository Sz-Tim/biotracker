/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package particle_track;

import java.util.Arrays;

/**
 * Only works for years 1993 - 2023
 * @author SA02MB
 */
public class ISO_datestr {
    private int day;
    private int month;
    private int year;
    private int[][] monthDays = new int[][]{{1,31},{2,28},{3,31},{4,30},{5,31},{6,30},{7,31},{8,31},{9,30},{10,31},{11,30},{12,31}};
    private final int[] leapYears = new int[]{1992,1996,2000,2004,2008,2012,2016,2020};
    
public ISO_datestr(int inDay, int inMonth, int inYear)
{
 this.day = inDay;
 this.month = inMonth;
 this.year  = inYear;
 
 for (int i = 0; i < this.leapYears.length; i++) {
            if(this.leapYears[i]==inYear)
            {monthDays[1][1]=29;}
        }
        
    
 }

public String getDateStr()
{
   String dateStr = String.format("%04d", this.year) + String.format("%02d", this.month) + String.format("%02d", this.day);   
   return dateStr;
}

public int getDateNum()
{
   // date integer with ref point year 1992, well why ever not?
   int yearDaysTot = (this.year - 1992)*365;
   int leapYearDaysTot = (int)Math.floor((this.year - 1992)/4) + 1;
   int monthDaysTot = 0;
   if (this.month >1){
   for(int i=0; i<=month-2;i++)
   {monthDaysTot = monthDaysTot + monthDays[i][1];       
   }}
   int dateNum =yearDaysTot + leapYearDaysTot + monthDaysTot + this.day;
   return dateNum;
}

public void addDay()
{
    if(this.day == this.monthDays[this.month -1][1])
       {this.day = 1;
        
        if(this.month==12)
        {this.month = 1;
         this.year++;
         boolean isLeapYear = false;
           for (int i = 0; i < this.leapYears.length; i++) {
            if(this.leapYears[i]==this.year)
            {isLeapYear=true;}
            }
           if(isLeapYear)
           {monthDays[1][1]=29;}
           else
           {monthDays[1][1]=28;}
        }           
        else
        {this.month++;}      
       }
    else{this.day++;} 
}

public void takeDay()
{
    if(this.day == 1)
       {if(this.month==1)
        {this.month = 12;
         this.day = this.monthDays[this.month -1][1];
         this.year--;
         boolean isLeapYear = false;
           for (int i = 0; i < this.leapYears.length; i++) {
            if(this.leapYears[i]==this.year)
            {isLeapYear=true;}
            }
           if(isLeapYear)
           {monthDays[1][1]=29;}
           else
           {monthDays[1][1]=28;}
        }           
        else
        {this.month--;
         this.day = this.monthDays[this.month][1];
        }      
       }
    else{this.day--;} 
}

}