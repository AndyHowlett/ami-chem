package org.xmlcml.ami2.chem;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;

import org.junit.Test;

public class ChemistryBuilderTest {

	@Test
	public void testDummy() {
		
	}
	
	@Test
	public void checkInChIs() throws IOException {
		BufferedReader r1 = new BufferedReader(new InputStreamReader(new FileInputStream("C:/workspace/inchisforevaluationtest.txt")));//trectrainingkekulised/inchisfinal.txt"
		String s = "";
		int i = 0;
		int numCorrect = 0;
		while ((s = r1.readLine()) != null) {
			/*if (i == 539 || i == 568 || i== 909) {
				i++;
			}*/
			try {
				//BufferedReader r2 = new BufferedReader(new InputStreamReader(new FileInputStream("C:/workspace/evaluationtestprimaryfragmentskekulised.xml.fromsvgsfromchembiodraw/molecule" + (i + 1) + ".svg.molecule0.txt")));//C:/workspace/evaluationtestprimaryfragmentskekulised.xml/molecule" + (i + 1) + ".svg.molecule0.txt")));//trectrainingkekulised.xml
				BufferedReader r2 = new BufferedReader(new InputStreamReader(new FileInputStream("C:/workspace/evaluationtestprimaryfragmentskekulised.xml/molecule" + (i + 1) + ".svg.molecule0.txt")));//trectrainingkekulised.xml
				boolean correct = r2.readLine().equals(s);
				//if (!correct) {
					System.out.println((i + 1) + " " + correct);
				//}
				numCorrect += (correct ? 1 : 0);
				if (!correct) {
					//FileUtils.copyFile(new File("C:/workspace/evaluationtestprimaryfragmentskekulised/molecule" + (i + 1) + ".svg"), new File("C:/workspace/evaluationtestprimaryfragmentskekulised.failuresforpixels/molecule" + (i + 1) + ".svg"));
				}
			} catch (Exception e) {
				//e.printStackTrace();
				System.out.println((i + 1) + " false");
				//FileUtils.copyFile(new File("C:/workspace/evaluationtestprimaryfragmentskekulised/molecule" + (i + 1) + ".svg"), new File("C:/workspace/evaluationtestprimaryfragmentskekulised.failuresforpixels/molecule" + (i + 1) + ".svg"));
			}
			i++;
		}
		System.out.println(numCorrect / ((float) i));
	}
	
}