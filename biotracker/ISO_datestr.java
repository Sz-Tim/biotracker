/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package biotracker;

/**
 * Only works for years 1993 - 2023 // updated SA04TS: leap years fixed; works until 2100
 *
 * @author SA02MB
 */
public class ISO_datestr {
    private int day;
    private int month;
    private int year;
    private final int[][] monthDays = new int[][]{{1, 31}, {2, 28}, {3, 31}, {4, 30}, {5, 31}, {6, 30}, {7, 31}, {8, 31}, {9, 30}, {10, 31}, {11, 30}, {12, 31}};

    public ISO_datestr(int inDay, int inMonth, int inYear) {
        this.day = inDay;
        this.month = inMonth;
        this.year = inYear;
        monthDays[1][1] = this.isLeapYear() ? 29 : 28;
    }

    public ISO_datestr(String dateString) {
        this.year = Integer.parseInt(dateString.substring(0, 4));
        this.month = Integer.parseInt(dateString.substring(4, 6));
        this.day = Integer.parseInt(dateString.substring(6, 8));
        monthDays[1][1] = this.isLeapYear() ? 29 : 28;
    }

    public ISO_datestr(ISO_datestr datestr) {
        this.year = datestr.getYear();
        this.month = datestr.getMonth();
        this.day = datestr.getDay();
        monthDays[1][1] = this.isLeapYear() ? 29 : 28;
    }

    public String getDateStr() {
        return String.format("%04d", this.year) + String.format("%02d", this.month) + String.format("%02d", this.day);
    }

    public int getDateNum() {
        // date integer with ref point year 1992, well why ever not?
        int yearDaysTot = (this.year - 1992) * 365;
        int leapYearDaysTot = (int) Math.floor((float) (this.year - 1992) / 4) + 1;
        int monthDaysTot = 0;
        if (this.month > 1) {
            for (int i = 0; i <= month - 2; i++) {
                monthDaysTot = monthDaysTot + monthDays[i][1];
            }
        }
        return yearDaysTot + leapYearDaysTot + monthDaysTot + this.day;
    }

    public void addDay() {
        if (this.day == this.monthDays[this.month - 1][1]) {
            this.day = 1;

            if (this.month == 12) {
                this.month = 1;
                this.year++;
                monthDays[1][1] = this.isLeapYear() ? 29 : 28;
            } else {
                this.month++;
            }
        } else {
            this.day++;
        }
    }

    public static ISO_datestr getTomorrow(ISO_datestr today) {
        ISO_datestr tomorrow = new ISO_datestr(today.getDay(), today.getMonth(), today.getYear());
        if (tomorrow.day == tomorrow.monthDays[tomorrow.month - 1][1]) {
            tomorrow.day = 1;

            if (tomorrow.month == 12) {
                tomorrow.month = 1;
                tomorrow.year++;
                tomorrow.monthDays[1][1] = tomorrow.isLeapYear() ? 29 : 28;
            } else {
                tomorrow.month++;
            }
        } else {
            tomorrow.day++;
        }
        return tomorrow;
    }

    public static ISO_datestr getYesterday(ISO_datestr today) {
        ISO_datestr yesterday = new ISO_datestr(today.getDay(), today.getMonth(), today.getYear());
        // If first of month, need to set day to be last day of previous month
        if (yesterday.day == 1) {
            if (yesterday.month == 1) {
                yesterday.month = 12;
                yesterday.year--;
                yesterday.monthDays[1][1] = yesterday.isLeapYear() ? 29 : 28;
            } else {
                yesterday.month--;
            }
            yesterday.day = yesterday.monthDays[yesterday.month - 1][1];
        } else {
            yesterday.day--;
        }
        return yesterday;
    }

    public void takeDay() {
        if (this.day == 1) {
            if (this.month == 1) {
                this.month = 12;
                this.day = this.monthDays[this.month - 1][1];
                this.year--;
                monthDays[1][1] = this.isLeapYear() ? 29 : 28;
            } else {
                this.month--;
                this.day = this.monthDays[this.month][1];
            }
        } else {
            this.day--;
        }
    }

    public int getYear() {
        return year;
    }

    public int getMonth() {
        return month;
    }

    public int getDay() {
        return day;
    }

    public boolean isLeapYear() {
        return this.year % 4 == 0;
    }

    public boolean isLaterThan(ISO_datestr dateCompare) {
        return this.getDateNum() > dateCompare.getDateNum();
    }

    /**
     * work out the date from an integer in format YYYYMMDD - Mike Bedington
     * method
     */
    public static int[] dateIntParse(int ymd) {
        double start_ymd_mod = ymd;
        int startYear = (int) Math.floor(start_ymd_mod / 10000);
        start_ymd_mod = start_ymd_mod - startYear * 10000;
        int startMonth = (int) Math.floor(start_ymd_mod / 100);
        start_ymd_mod = start_ymd_mod - startMonth * 100;
        int startDay = (int) start_ymd_mod;

        return new int[]{startDay, startMonth, startYear};
    }
    public static int[] dateIntParse(ISO_datestr datestr) {
        return new int[]{datestr.getDay(), datestr.getMonth(), datestr.getYear()};
    }

    @Override
    public String toString() {
        return this.year + "-" + (this.month < 10 ? "0" : "") + this.month + "-" + (this.day < 10 ? "0" : "") + this.day;
    }
    public String toStringShort() {
        return this.year + (this.month < 10 ? "0" : "") + this.month + (this.day < 10 ? "0" : "") + this.day;
    }
}