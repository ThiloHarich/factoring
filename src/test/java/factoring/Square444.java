package factoring;

import org.junit.Test;

/**
 * Created by Thilo Harich on 01.10.2017.
 */
public class Square444 {

    @Test
    public void chek444 ()
    {
        for(int i = 1;i<10000;i++)
        {
            if ((i * i) % 10000 == 4444)
                System.out.println(i);
        }
    }
}
