import java.io.*;
import java.util.zip.*;
import java.util.*;

public class MatrixReader{
    public static double [][] readMatrix(String file){
        try{
            FileInputStream fileInStream = new FileInputStream(new File(file));
            GZIPInputStream stream = new GZIPInputStream(fileInStream );
            BufferedReader reader = new BufferedReader(new InputStreamReader(stream), 1024*1024*16);
            List<double[]> data = new ArrayList<double[]>();
            int len = 0;
            while (reader.ready()){
                String line = reader.readLine();
                String tokens[] = line.split(",");
                if (len != tokens.length && len > 0)
                    System.out.println(line + " " + len);
                len = tokens.length;
                double []a = new double[tokens.length];
                for (int i=0; i<a.length; i++)
                    a[i] = Double.parseDouble(tokens[i]);
                data.add(a);
            }
            double [][] res = new double[data.size()][];
            for (int i=0; i<res.length; i++)
                res[i] = data.get(i);
            return res;
        } catch(IOException ex){
            ex.printStackTrace();
            return null;
        }

    }

    public static void main(String[] args){
        readMatrix(args[0]);
    }
}

