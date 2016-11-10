package org.xmlcml.ami2.chem;

import java.awt.Graphics2D;
import java.awt.font.GlyphVector;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.concurrent.TimeoutException;

import javax.imageio.ImageIO;

import nu.xom.Nodes;

import org.apache.commons.math.complex.Complex;
import org.apache.log4j.Logger;
import org.xmlcml.ami2.chem.Joinable.JoinPoint;
import org.xmlcml.ami2.chem.JoinableText.AreInSameStringDetector;
import org.xmlcml.ami2.chem.svg.SVGContainerNew;
import org.xmlcml.diagrams.DiagramAnalyzer;
import org.xmlcml.diagrams.OCRManager;
import org.xmlcml.euclid.Angle;
import org.xmlcml.euclid.Angle.Units;
import org.xmlcml.euclid.Line2;
import org.xmlcml.euclid.Real;
import org.xmlcml.euclid.Real2;
import org.xmlcml.euclid.Real2Array;
import org.xmlcml.euclid.Real2Range;
import org.xmlcml.euclid.RealRange;
import org.xmlcml.euclid.Transform2;
import org.xmlcml.euclid.Vector2;
import org.xmlcml.graphics.svg.SVGCircle;
import org.xmlcml.graphics.svg.SVGConstants;
import org.xmlcml.graphics.svg.SVGElement;
import org.xmlcml.graphics.svg.SVGG;
import org.xmlcml.graphics.svg.SVGImage;
import org.xmlcml.graphics.svg.SVGLine;
import org.xmlcml.graphics.svg.SVGPath;
import org.xmlcml.graphics.svg.SVGPolygon;
import org.xmlcml.graphics.svg.SVGSVG;
import org.xmlcml.graphics.svg.SVGTSpan;
import org.xmlcml.graphics.svg.SVGText;
import org.xmlcml.svgbuilder.geom.SimpleBuilder;

import com.google.common.collect.ImmutableCollection;
import com.google.common.collect.UnionFind;
import com.google.common.util.concurrent.UncheckedTimeoutException;
//import net.sourceforge.tess4j.Tesseract;
//import net.sourceforge.tess4j.TesseractException;

/**
 * Builds higher-level primitives from SVGPaths, SVGLines, etc. to create SVG objects 
 * such as TramLine and (later) Arrow.
 * 
 * <p>SimpleBuilder's main function is to:
 * <ul>
 * <li>Read a raw SVG object and make lists of SVGPath and SVGText (and possibly higher levels ones
 * if present).</li>
 * <li>Turn SVGPaths into SVGLines , etc..</li>
 * <li>Identify Junctions (line-line, line-text, and probably more).</li>
 * <li>Join lines where they meet into higher level objects (TramLines, SVGRect, crosses, arrows, etc.).</li>
 * <li>Create topologies (e.g. connection of lines and Junctions).</li>
 * </ul>
 * 
 * SimpleBuilder uses the services of the org.xmlcml.graphics.svg.path package and may later use
 * org.xmlcml.graphics.svg.symbol.
 * </p>
 * 
 * <p>Input may either be explicit SVG primitives (e.g. &lt;svg:rect&gt;, &lt;svg:line&gt;) or 
 * implicit ones (&lt;svg:path&gt;) that can be interpreted as the above. The input may be either or
 * both - we can't control it. The implicit get converted to explicit and then merged with the 
 * explicit:
 * <pre>
 *    paths-> implicitLineList + rawLinelist -> explicitLineList 
 * </pre>
 * </p>
 * 
 * <h3>Strategy</h3>
 * <p>createHigherLevelPrimitives() carries out the complete chain from svgRoot to the final
 * primitives. Each step tests to see whether the result of the previous is null.
 * If so it creates a non-null list and fills it if possible.</p>
 * 
 * @author pm286
 */
public class ChemistryBuilder extends SimpleBuilder {
	
	private final static Logger LOG = Logger.getLogger(ChemistryBuilder.class);
	
	private ChemistryBuilderParameters parameters = new ChemistryBuilderParameters();

	protected HigherPrimitives higherPrimitives;

	protected SVGContainerNew input;

	private List<JoinableText> atomLabelTexts;
	private Map<Real2Range, Integer> atomLabelPositionsAndNumbers;
	
	List<WedgeBond> wedgeBonds = new ArrayList<WedgeBond>();
	
	private double scale = 1;

	AmbiguityHandler ambiguityHandler = new AmbiguityHandler(this);

	public ChemistryBuilder(SVGContainerNew svgRoot, long timeout, ChemistryBuilderParameters parameters) {
		super((SVGElement) svgRoot.getElement(), timeout);
		input = svgRoot;
		this.parameters = parameters;
	}

	public ChemistryBuilder(SVGContainerNew svgRoot, ChemistryBuilderParameters parameters) {
		super((SVGElement) svgRoot.getElement());
		input = svgRoot;
		this.parameters = parameters;
	}
	
	public ChemistryBuilder(SVGElement svgRoot, long timeout, ChemistryBuilderParameters parameters) {
		super(svgRoot, timeout);
		this.parameters = parameters;
	}

	public ChemistryBuilder(SVGElement svgRoot, ChemistryBuilderParameters parameters) {
		super(svgRoot);
		this.parameters = parameters;
	}
	
	public ChemistryBuilder(SVGContainerNew svgRoot, long timeout) {
		super((SVGElement) svgRoot.getElement(), timeout);
		input = svgRoot;
		parameters = new ChemistryBuilderParameters();
	}

	public ChemistryBuilder(SVGContainerNew svgRoot) {
		super((SVGElement) svgRoot.getElement());
		input = svgRoot;
		parameters = new ChemistryBuilderParameters();
	}
	
	public ChemistryBuilder(SVGElement svgRoot, long timeout) {
		super(svgRoot, timeout);
		parameters = new ChemistryBuilderParameters();
	}

	public ChemistryBuilder(SVGElement svgRoot) {
		super(svgRoot);
		parameters = new ChemistryBuilderParameters();
	}
	
	public ChemistryBuilderParameters getParameters() {
		return parameters;
	}
	
	void setParameters(ChemistryBuilderParameters parameters) {
		this.parameters = parameters;
	}

	public SVGContainerNew getInputContainer() {
		return input;
	}
	
	/**
	 * Complete processing chain for low-level SVG into high-level SVG and non-SVG primitives such as double bonds.
	 * <p>
	 * Creates junctions.
	 * <p>
	 * Runs createDerivedPrimitives().
	 * 
	 * @throws TimeoutException 
	 */
	public void createHigherPrimitives() {
		if (higherPrimitives == null) {
			startTiming();
			//replaceTextVectorsWithText();
			createDerivedPrimitives();
			replaceTextImagesWithText();
			splitMultiCharacterTexts();
			higherPrimitives = new HigherPrimitives();
			higherPrimitives.addSingleLines(derivedPrimitives.getLineList());
			handleShortLines();
			createUnsaturatedBondLists();
			//createWords();
			createJunctions();
		}
	}

	@Override
	protected void removeNearDuplicateAndObscuredPrimitives() {
		double scale = parameters.setStandardBondLengthFromSVG(derivedPrimitives.getLineList());
		nearDuplicateLineRemovalDistance *= (scale / this.scale);
		nearDuplicatePolygonRemovalDistance *= (scale / this.scale);
		minimumCutObjectGap *= (scale / this.scale);
		maximumCutObjectGap *= (scale / this.scale);
		this.scale = scale;
		super.removeNearDuplicateAndObscuredPrimitives();
	}
	
	private void splitMultiCharacterTexts() {
		Iterator<SVGText> it = derivedPrimitives.getTextList().iterator();
		List<SVGText> newTexts = new ArrayList<SVGText>();
		while (it.hasNext()) {
			SVGText text = it.next();
			String string = text.getText();
			List<SVGTSpan> spanList = new ArrayList<SVGTSpan>();
			double totalWidth = 0;
			if (string == null) {
				Nodes spans = text.query("svg:tspan", SVGSVG.SVG_XPATH);
				for (int i = 0; i < spans.size(); i++) {
					SVGTSpan span = (SVGTSpan) spans.get(i);
					spanList.add(span);
					GlyphVector v = span.getGlyphVector();
					totalWidth += v.getLogicalBounds().getWidth();
					if (span.getAttributeValue("dx") != null) {
						totalWidth += Double.parseDouble(span.getAttributeValue("dx"));
					}
				}
			} else {
				spanList.add(new SVGTSpan(text));
				totalWidth = text.getGlyphVector().getLogicalBounds().getWidth();
			}
			it.remove();
			double previousX = text.getX() - ("end".equals(text.getAttributeValue("text-anchor")) ? totalWidth : 0);
			double previousY = text.getY();
			for (SVGTSpan span : spanList) {
				if (span.getX() != 0.0) {
					previousX = span.getX() - ("end".equals(text.getAttributeValue("text-anchor")) ? totalWidth : 0);
				}
				if (span.getY() != 0.0) {
					previousY = span.getY();
				}
				GlyphVector glyphVector = span.getGlyphVector();
				String spanText = span.getText();
				for (int i = 0; i < spanText.length(); i++) {
					String substring = spanText.substring(i, i + 1);
					SVGText newText = new SVGText(new Real2(0, 0), substring);
					newTexts.add(newText);
					newText.copyAttributesFrom(span);
					double dX = 0;
					if (span.getAttributeValue("dx") != null) {
						dX = Double.parseDouble(span.getAttributeValue("dx"));
						newText.removeAttribute(newText.getAttribute("dx"));
					}
					double dY = 0;
					if (span.getAttributeValue("dy") != null) {
						dY = Double.parseDouble(span.getAttributeValue("dy"));
						newText.removeAttribute(newText.getAttribute("dy"));
					}
					newText.setX(previousX + dX + glyphVector.getGlyphPosition(i).getX());
					newText.setY(previousY + dY + glyphVector.getGlyphPosition(i).getY());
				}
				previousX += glyphVector.getGlyphPosition(glyphVector.getNumGlyphs()).getX();
				if (span.getAttributeValue("dy") != null) {
					previousY += Double.parseDouble(span.getAttributeValue("dy"));
				}
			}
		}
		derivedPrimitives.getTextList().addAll(newTexts);
	}

	public BufferedImage flipHorizontally(BufferedImage img) {
		int w = img.getWidth();
		int h = img.getHeight();
		BufferedImage dimg = new BufferedImage(w, h, img.getColorModel().getTransparency());
		Graphics2D g = dimg.createGraphics();
		g.drawImage(img, 0, 0, w, h, 0, h, w, 0, null);
		g.dispose();
		return dimg;
	}

	private File getImageFileFromSVGImage(SVGImage image) {
		String filename = image.getAttributeValue("href", SVGConstants.XLINK_NS);
		File testFile = new File(filename);
		return (testFile.isAbsolute() ? testFile : new File(input.getFile().getParentFile().getAbsolutePath() + "/" + filename));
	}
	
	private void replaceTextVectorsWithText() {
		List<SVGPath> paths = getRawPrimitives().getPathList();
		Real2Range range = new Real2Range();
		for (SVGPath path : paths) {
			range.plusEquals(path.getBoundingBox());
		}
		BufferedImage im = new BufferedImage((int) range.getXRange().getRange(), (int) range.getYRange().getRange(), BufferedImage.TYPE_INT_ARGB);
		Transform2 translate = new Transform2(new Vector2(-range.getXMin(), -range.getYMin()));
		for (SVGPath path : paths) {
			path.applyTransform(translate);
			path.draw(im.createGraphics());
		}
		DiagramAnalyzer analyzer = new DiagramAnalyzer();
		analyzer.setImage(im);
		SVGSVG out = analyzer.convertPixelsToSVG();
	}
	
	private void replaceTextImagesWithText() {
		if (rawPrimitives.getImageList().size() == 0) {
			return;
		}
		Set<Complex> done = new HashSet<Complex>();
		OCRManager manager = new OCRManager();
		
		for (SVGImage image : rawPrimitives.getImageList()) {
			try {
				checkTime("Took too long to convert images to text");
				image.applyTransformAttributeAndRemove();
				if (image.getWidth() > parameters.getMaximumImageElementWidthForOCR()) {
					continue;
				}
				File file = getImageFileFromSVGImage(image);
				BufferedImage bufferedImage = flipHorizontally(ImageIO.read(file));
				/*Tesseract tess = Tesseract.getInstance();
				try {
					s = tess.doOCR(im);
				} catch (TesseractException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}*/
				SVGText text = manager.scan(bufferedImage, new Real2Range(new RealRange(image.getX(), image.getX() + image.getWidth()), new RealRange(Math.min(image.getY(), image.getY() + image.getHeight()), Math.max(image.getY(), image.getY() + image.getHeight()))), parameters.getBlackThreshold(), parameters.getMaximumOCRError());
				if (text != null) {
					image.getParent().replaceChild(image, text);
					//text.copyAttributesFrom(image);
					derivedPrimitives.getTextList().add(text);
					derivedPrimitives.getImageList().remove(image);
				}
				
				if (!done.add(new Complex(image.getX(), image.getY())) || bufferedImage.getWidth() < parameters.getMimimumImageWidthForOCR()) {
					derivedPrimitives.getImageList().remove(image);
					continue;
				}
			} catch (IOException e) {
				System.err.println("Error handling image within SVG file - it's probably embedded in base 64, but it should be linked to and stored separately");
				e.printStackTrace();
			} catch (Exception e) {
				//TODO handle other images
			}

		}
		
		manager.handleAmbiguousTexts(parameters.getTextCoordinateTolerance(), parameters.getAllowedFontSizeVariation());
	}
	
	/*private void convertImagesOfTextToText() {
		for (SVGImage image : rawPrimitives.getImageList()) {
			String path = image.getAttributeValue("href", "xlink");
		}
	}*/

	/*private void createWords() {
		TextStructurer t = new TextStructurer(derivedPrimitives.getTextList());
		//List<RawWords> lines = t.createRawWordsList();
		List<ScriptLine> lines = t.getScriptedLineList();
		//j.get(0).getScriptWordList().get(0).createSuscriptTextLineList();
		List<JoinableScriptWord> words = new ArrayList<JoinableScriptWord>();
		higherPrimitives.setWordsList(words);
		for (ScriptLine line : lines) {
			for (ScriptWord word : line.getScriptWordList()) {
				words.add(new JoinableScriptWord(word));
			}
		}
	}*/

	private void handleShortLines() {
		List<HatchedBond> hatchList = new ArrayList<HatchedBond>();
		higherPrimitives.setHatchedBondList(hatchList);
		List<SVGLine> smallLines = new ArrayList<SVGLine>();
		for (SVGLine l : derivedPrimitives.getLineList()) {
			if (l.getXY(0).getDistance(l.getXY(1)) < parameters.getHatchLineMaximumLength() && l.getXY(0).getDistance(l.getXY(1)) > 0) {//TODO l.getLength() < hatchLineMaximumLength) {
				smallLines.add(l);
			}
		}
		higherPrimitives.setChargeList(new ArrayList<Charge>());
		if (smallLines.size() == 0) {
			ambiguityHandler.mutuallyExclusiveShortLineTriples = new ArrayList<AmbiguityHandler.MutuallyExclusiveShortLineTriple>();
			ambiguityHandler.mutuallyExclusiveShortLinePairTriples = new ArrayList<AmbiguityHandler.MutuallyExclusiveShortLinePairTriple>();
			ambiguityHandler.mutuallyExclusiveShortLineTripleTriples = new ArrayList<AmbiguityHandler.MutuallyExclusiveShortLineTripleTriple>();
			return;
		}
		
		sewSmallLines(smallLines);
		
		UnionFind<SVGLine> hatchedBonds = UnionFind.create(smallLines);
		for (int i = 0; i < smallLines.size(); i++) {
			SVGLine firstLine = smallLines.get(i);
			for (int j = i + 1; j < smallLines.size(); j++) {
				checkTime("Took too long to handle short lines");
				SVGLine secondLine = smallLines.get(j);
				Double dist = firstLine.calculateUnsignedDistanceBetweenLines(secondLine, new Angle((firstLine.getLength() < parameters.getTinyHatchLineMaximumLength() || secondLine.getLength() < parameters.getTinyHatchLineMaximumLength() ? parameters.getMaximumAngleForParallelIfOneLineIsTiny() : parameters.getMaximumAngleForParallel()), Units.RADIANS));
				if (dist != null && dist < parameters.getHatchLinesMaximumSpacing() && dist > parameters.getHatchLinesMinimumSpacing() && (firstLine.overlapsWithLine(secondLine, parameters.getLineOverlapEpsilon()) || secondLine.overlapsWithLine(firstLine, parameters.getLineOverlapEpsilon()))) {
					try {
						hatchedBonds.union(firstLine, secondLine);
					} catch (IllegalArgumentException e) {
						
					}
				}
				if ((firstLine.isHorizontal(parameters.getFlatLineEpsilon()) || secondLine.isHorizontal(parameters.getFlatLineEpsilon())) && firstLine.overlapsWithLine(secondLine, parameters.getLineOverlapEpsilon()) && secondLine.overlapsWithLine(firstLine, parameters.getLineOverlapEpsilon()) && firstLine.getEuclidLine().isPerpendicularTo(secondLine.getEuclidLine(), new Angle(parameters.getPlusChargeAngleTolerance(), Units.DEGREES))) {
					hatchedBonds.remove(firstLine);
					hatchedBonds.remove(secondLine);
					higherPrimitives.getLineList().remove(firstLine);
					higherPrimitives.getLineList().remove(secondLine);
					higherPrimitives.getLineChargeList().add(new Charge(parameters, firstLine, secondLine));
				}
			}
		}
		handleShortLines(hatchedBonds);
	}

	private void sewSmallLines(List<SVGLine> smallLines) {
		Iterator<SVGLine> it = smallLines.iterator();
		int i = 0;
		outer: for (SVGLine line1 = it.next(); it.hasNext(); line1 = it.next()) {
			for (int j = i + 1; j < smallLines.size(); j++) {
				SVGLine line2 = smallLines.get(j);
				boolean test1 = line1.getXY(0).isEqualTo(line2.getXY(0), 1.5);
				boolean test2 = line1.getXY(0).isEqualTo(line2.getXY(1), 1.5);
				boolean test3 = line1.getXY(1).isEqualTo(line2.getXY(0), 1.5);
				boolean test4 = line1.getXY(1).isEqualTo(line2.getXY(1), 1.5);
				if (test1 || test2 || test3 || test4) {
					if (line1.isParallelOrAntiParallelTo(line2, new Angle(parameters.getMaximumAngleForParallelIfOneLineIsTiny(), Units.RADIANS))) {
						Real2 common = line1.getCommonEndPoint(line2, 1.5);
						SVGLine newLine = new SVGLine(line2.getOtherPoint(common, 1.5), line1.getOtherPoint(common, 1.5));
						if (newLine.isParallelOrAntiParallelTo(line1, new Angle(parameters.getMaximumAngleForParallel(), Units.RADIANS)) && newLine.isParallelOrAntiParallelTo(line2, new Angle(parameters.getMaximumAngleForParallel(), Units.RADIANS))) {
							it.remove();
							higherPrimitives.getLineList().remove(line1);
							line2.setXY(newLine.getXY(0), 0);
							line2.setXY(newLine.getXY(1), 1);
							continue outer;
						}
					}
				}
			}
			i++;
		}
	}

	private void handleShortLines(UnionFind<SVGLine> disjointSets) {
		final double threshold = parameters.getThresholdForOrderingCheckForHatchedBonds();
		ambiguityHandler.mutuallyExclusiveShortLineTriples = new ArrayList<AmbiguityHandler.MutuallyExclusiveShortLineTriple>();
		ambiguityHandler.mutuallyExclusiveShortLinePairTriples = new ArrayList<AmbiguityHandler.MutuallyExclusiveShortLinePairTriple>();
		ambiguityHandler.mutuallyExclusiveShortLineTripleTriples = new ArrayList<AmbiguityHandler.MutuallyExclusiveShortLineTripleTriple>();
		List<HatchedBond> hatchList = higherPrimitives.getHatchedBondList();
		set: for (Set<SVGLine> set : disjointSets.snapshot()) {
			ArrayList<SVGLine> lines1 = new ArrayList<SVGLine>(set);
			ArrayList<SVGLine> lines2 = new ArrayList<SVGLine>(set);
			Collections.sort(lines1, new Comparator<SVGLine>(){
				public int compare(SVGLine i, SVGLine j) {
					Real2 firstOfI = (i.getXY(0).getX() < i.getXY(1).getX() ? i.getXY(0) : i.getXY(1));
					Real2 firstOfJ = (j.getXY(0).getX() < j.getXY(1).getX() ? j.getXY(0) : j.getXY(1));
					return (Real.isEqual(firstOfI.getX(), firstOfJ.getX(), threshold) ? Double.compare(firstOfI.getY(), firstOfJ.getY()) : Double.compare(firstOfI.getX(), firstOfJ.getX()));
				}});
			Collections.sort(lines2, new Comparator<SVGLine>(){
				public int compare(SVGLine i, SVGLine j) {
					Real2 firstOfI = (i.getXY(0).getY() < i.getXY(1).getY() ? i.getXY(0) : i.getXY(1));
					Real2 firstOfJ = (j.getXY(0).getY() < j.getXY(1).getY() ? j.getXY(0) : j.getXY(1));
					return (Real.isEqual(firstOfI.getY(), firstOfJ.getY(), threshold) ? Double.compare(firstOfI.getX(), firstOfJ.getX()) : Double.compare(firstOfI.getY(), firstOfJ.getY()));
				}});
			ArrayList<SVGLine> lines3 = (ArrayList<SVGLine>) lines1.clone();
			Collections.reverse(lines3);
			ArrayList<SVGLine> lines4 = (ArrayList<SVGLine>) lines1.clone();
			ArrayList<SVGLine> lines5 = (ArrayList<SVGLine>) lines1.clone();
			Collections.sort(lines4, new Comparator<SVGLine>(){
				public int compare(SVGLine i, SVGLine j) {
					Real2 secondOfI = (i.getXY(0).getX() >= i.getXY(1).getX() ? i.getXY(0) : i.getXY(1));
					Real2 secondOfJ = (j.getXY(0).getX() >= j.getXY(1).getX() ? j.getXY(0) : j.getXY(1));
					return (Real.isEqual(secondOfI.getX(), secondOfJ.getX(), threshold) ? Double.compare(secondOfI.getY(), secondOfJ.getY()) : Double.compare(secondOfI.getX(), secondOfJ.getX()));
				}});
			Collections.sort(lines5, new Comparator<SVGLine>(){
				public int compare(SVGLine i, SVGLine j) {
					Real2 secondOfI = (i.getXY(0).getY() >= i.getXY(1).getY() ? i.getXY(0) : i.getXY(1));
					Real2 secondOfJ = (j.getXY(0).getY() >= j.getXY(1).getY() ? j.getXY(0) : j.getXY(1));
					return (Real.isEqual(secondOfI.getY(), secondOfJ.getY(), threshold) ? Double.compare(secondOfI.getX(), secondOfJ.getX()) : Double.compare(secondOfI.getY(), secondOfJ.getY()));
				}});
			ArrayList<SVGLine> lines6 = (ArrayList<SVGLine>) lines4.clone();
			Collections.reverse(lines6);
			ArrayList<SVGLine> lines7 = (ArrayList<SVGLine>) lines1.clone();
			ArrayList<SVGLine> lines8 = (ArrayList<SVGLine>) lines1.clone();
			Collections.sort(lines7, new Comparator<SVGLine>(){
				public int compare(SVGLine i, SVGLine j) {
					return (Real.isEqual(i.getMidPoint().getX(), j.getMidPoint().getX(), threshold) ? Double.compare(i.getMidPoint().getY(), j.getMidPoint().getY()) : Double.compare(i.getMidPoint().getX(), j.getMidPoint().getX()));
				}});
			Collections.sort(lines8, new Comparator<SVGLine>(){
				public int compare(SVGLine i, SVGLine j) {
					return (Real.isEqual(i.getMidPoint().getY(), j.getMidPoint().getY(), threshold) ? Double.compare(i.getMidPoint().getX(), j.getMidPoint().getX()) : Double.compare(i.getMidPoint().getY(), j.getMidPoint().getY()));
				}});
			ArrayList<SVGLine> lines9 = (ArrayList<SVGLine>) lines7.clone();
			Collections.reverse(lines9);
			boolean firstEndPointsAndSecondEndPointsOrdered = ((lines1.equals(lines2) || lines3.equals(lines2)) && (lines4.equals(lines5) || lines6.equals(lines5)));
			boolean firstEndPointsAndMidPointsOrdered = ((lines1.equals(lines2) || lines3.equals(lines2)) && (lines7.equals(lines8) || lines9.equals(lines8)));
			boolean secondEndPointsAndMidPointsOrdered = ((lines4.equals(lines5) || lines6.equals(lines5)) && (lines7.equals(lines8) || lines9.equals(lines8)));
			if (firstEndPointsAndSecondEndPointsOrdered || firstEndPointsAndMidPointsOrdered || secondEndPointsAndMidPointsOrdered) {
				ArrayList<SVGLine> lines;
				if (firstEndPointsAndSecondEndPointsOrdered || firstEndPointsAndMidPointsOrdered) {
					lines = lines1;
				} else {
					lines = lines4;
				}
				try {
					double lengthFirst = lines.get(0).getLength();
					double lengthLast = lines.get(lines.size() - 1).getLength();
					double shortest = Double.MAX_VALUE;
					double longest = Double.MIN_VALUE;
					boolean allLongerFirst = true;
					boolean allShorterFirst = true;
					boolean allLongerLast = true;
					boolean allShorterLast = true;
					for (int i = 0; i < lines.size(); i++) {
						double length = lines.get(i).getLength();
						if (length < shortest) {
							shortest = length;
						}
						if (length > longest) {
							longest = length;
						}
						allLongerFirst &= length >= lengthFirst;
						allShorterFirst &= length <= lengthFirst;
						allLongerLast &= length >= lengthLast;
						allShorterLast &= length <= lengthLast;
					}
					if (longest - shortest > parameters.getLengthTolerance() && !allLongerFirst && !allShorterFirst && !allLongerLast && !allShorterLast) {
						continue set;
					}
				} catch (IndexOutOfBoundsException e) {
					
				}
				HatchedBond hatchedBond = new HatchedBond(parameters, lines);
				hatchList.add(hatchedBond);
				if (lines.size() > 3) {
					higherPrimitives.getLineList().removeAll(lines);
				} else if (lines.size() == 3) {
					ambiguityHandler.mutuallyExclusiveShortLineTripleTriples.add(new AmbiguityHandler.MutuallyExclusiveShortLineTripleTriple(hatchedBond, lines.get(0), lines.get(1), lines.get(2)));
				} else if (lines.size() == 2) {
					ambiguityHandler.mutuallyExclusiveShortLinePairTriples.add(new AmbiguityHandler.MutuallyExclusiveShortLinePairTriple(hatchedBond, lines.get(0), lines.get(1)));
				} else {
					Charge charge = null;
					if (lines.get(0).isHorizontal(parameters.getFlatLineEpsilon()) && !lines.get(0).isVertical(parameters.getFlatLineEpsilon())) {
						charge = new Charge(parameters, lines);
						higherPrimitives.getLineChargeList().add(charge);
					}
					ambiguityHandler.mutuallyExclusiveShortLineTriples.add(new AmbiguityHandler.MutuallyExclusiveShortLineTriple(hatchedBond, charge, lines.get(0)));
				}
			}
		}
	}
	
	private void createJunctions() {
		createJoinableList();
		List<Joinable> joinables = higherPrimitives.getJoinableList();

		List<JoinPoint> joinPoints = extractAtomLabelsAndGetRemainingJoinPoints(joinables);
		
		UnionFind<JoinPoint> joinPointsGroupedIntoJunctions = UnionFind.create(joinPoints);
		//int deleted = 0;
		attemptToJoinListOfJoinables(joinables, joinPointsGroupedIntoJunctions);
		
		try {
			ambiguityHandler.handleAmbiguities(joinPointsGroupedIntoJunctions);
		} catch (IllegalArgumentException e) {
			joinPoints.clear();
			LOG.debug("Processing failed as the diagram was too complex");
		}
		
		List<Junction> junctions = new ArrayList<Junction>();
			
		if (joinPoints.size() != 0) {
			for (Set<JoinPoint> junctionJoinPoints : joinPointsGroupedIntoJunctions.snapshot()) {
				//if (junctionJoinPoints.size() != 1) {
					/*Set<Joinable> junctionJoinables = new HashSet<Joinable>();
					Set<JoinPoint> newJunctionJoinPoints = new HashSet <JoinPoint>();
					for (JoinPoint point : junctionJoinPoints) {
						if (junctionJoinables.add(point.getJoinable())) {
							newJunctionJoinPoints.add(point);
						} else {
							newJunctionJoinPoints.removeAll(point.getJoinable().getJoinPoints());
							removeJoinable(point.getJoinable());
							//junctionJoinables.remove(point.getJoinable());
						}
					}
					for (Joinable j1 : junctionJoinables) {
						int numberParallel = 0;
						for (Joinable j2 : junctionJoinables) {
							if (Joinable.areParallel(j1, j2)) {
								numberParallel++;
								if (numberParallel == 3) {
									for (JoinPoint p : newJunctionJoinPoints) {
										junctions.add(new Junction(p));
									}
									continue junction;
								}
							}
						}
					}
					junctions.add(new Junction(newJunctionJoinPoints));*/
				//}
				junctions.add(new Junction(junctionJoinPoints));
			}
		}
		higherPrimitives.setJunctionList(junctions);
				
				/*JoinPoint commonPoint = joinablei.getIntersectionPoint(joinablej);
				if (commonPoint != null) {
					Junction junction = new Junction(joinablei, joinablej, commonPoint);
					rawJunctionList.add(junction);
					String junctAttVal = "junct"+"."+rawJunctionList.size();
					junction.addAttribute(new Attribute(SVGElement.ID, junctAttVal));
					if (junction.getCoordinates() == null && commonPoint.getPoint() != null) {
						junction.setCoordinates(commonPoint.getPoint());
					}
					LOG.debug("junct: "+junction.getId()+" between " + joinablei.getClass() + " and " + joinablej.getClass() + " with coords "+junction.getCoordinates()+" "+commonPoint.getPoint());
				}
			}
		}*/
		
		/*createRawJunctionList();
		List<Junction> junctionList = new ArrayList<Junction>(higherPrimitives.getRawJunctionList());
		for (int i = junctionList.size() - 1; i > 0; i--) {
			Junction labile = junctionList.get(i);
			for (int j = 0; j < i; j++) {
				Junction fixed = junctionList.get(j);
				if (fixed.containsCommonPoints(labile)) {
					labile.transferDetailsTo(fixed);
					junctionList.remove(i);
					break;
				}
			}
		}
		higherPrimitives.setMergedJunctionList(junctionList);*/
	}

	void attemptToJoinListOfJoinables(List<Joinable> joinables, UnionFind<JoinPoint> joinPointsGroupedIntoJunctions) {
		List<JoinableText> texts = new ArrayList<JoinableText>();
		for (int i = 0; i < joinables.size() - 1; i++) {
			Joinable joinableI = joinables.get(i);
			if (!(joinableI instanceof JoinableText)) {
				for (int j = i + 1; j < joinables.size(); j++) {
					Joinable joinableJ = joinables.get(j);
					if (!(joinableJ instanceof JoinableText)) {
						checkTime("Took too long to determine what is joined to what");
						
						//System.out.println(joinableI + "\n" + joinableJ);
						
						joinPointsGroupedIntoJunctions.unionAll(getListOfOverlappingJoinPointsForJoinables(joinPointsGroupedIntoJunctions, joinableI, joinableJ));
					}
				}
			} else {
				texts.add((JoinableText) joinableI);
			}
		}
		if (joinables.size() > 0) {
			if (joinables.get(joinables.size() - 1) instanceof JoinableText) {
				texts.add((JoinableText) joinables.get(joinables.size() - 1));
			}
			attemptToJoinTexts(texts, joinables, joinPointsGroupedIntoJunctions);
		}
	}

	private void attemptToJoinTexts(List<JoinableText> texts, List<Joinable> joinables, UnionFind<JoinPoint> joinPointsGroupedIntoJunctions) {
		for (int i = 0; i < texts.size() - 1; i++) {
			JoinableText textI = texts.get(i);
			for (int j = i + 1; j < texts.size(); j++) {
				JoinableText textJ = texts.get(j);
				checkTime("Took too long to determine what is joined to what");
				if (!getListOfOverlappingJoinPointsForJoinables(joinPointsGroupedIntoJunctions, textI, textJ).isEmpty()) {
					System.out.println(textI.getSVGElement().getText() + " " + textJ.getSVGElement().getText());					
				}
				joinPointsGroupedIntoJunctions.unionAll(getListOfOverlappingJoinPointsForJoinables(joinPointsGroupedIntoJunctions, textI, textJ));
			}
		}
		for (JoinableText text : texts) {
			Set<JoinPoint> joinPoints = new HashSet<JoinPoint>();
			for (Joinable joinable : joinables) {
				if (!(joinable instanceof JoinableText)) {
					checkTime("Took too long to determine what is joined to what");
					List<JoinPoint> overlap = getListOfOverlappingJoinPointsForJoinables(joinPointsGroupedIntoJunctions, text, joinable);
					joinPoints.addAll(overlap);
					for (JoinPoint j : overlap) {
						if (!(j.getJoinable() instanceof JoinableText)) {
							joinPoints.addAll(joinPointsGroupedIntoJunctions.getObjectsInPartitionOf(j));
						}
					}
				}
			}
			List<JoinPoint> actualJoinPoints = new ArrayList<JoinPoint>();
			i: for (JoinPoint joinPointI : joinPoints) {
				for (JoinPoint joinPointJ : joinPoints) {
					if (joinPointI != joinPointJ && joinPointsGroupedIntoJunctions.getObjectsInPartitionOf(joinPointI).contains(joinPointJ)) {
						//continue i;
					}
				}
				actualJoinPoints.add(joinPointI);
			}
			joinPointsGroupedIntoJunctions.unionAll(actualJoinPoints);
		}
	}

	private List<JoinPoint> getListOfOverlappingJoinPointsForJoinables(UnionFind<JoinPoint> joinPointsGroupedIntoJunctions, Joinable joinableI, Joinable joinableJ) {
		Set<JoinPoint> overlapSet = joinableI.overlapWith(joinableJ);
		if (overlapSet != null) {
			List<JoinPoint> overlapList = new ArrayList<JoinPoint>(overlapSet);
			if (!joinPointsGroupedIntoJunctions.contains(overlapList.get(0)) || !joinPointsGroupedIntoJunctions.contains(overlapList.get(1)) || (overlapList.size() > 2 && !joinPointsGroupedIntoJunctions.contains(overlapList.get(2))) || (overlapList.size() > 3 && !joinPointsGroupedIntoJunctions.contains(overlapList.get(3)))) {
				return new ArrayList<JoinPoint>();
			}
			if (ambiguityHandler.mutuallyExclusive(joinableI, joinableJ)) {
				return new ArrayList<JoinPoint>();
			}
			if (joinableI instanceof JoinableText && joinableJ instanceof JoinableText) {
				if (JoinableText.doTextsJoin((JoinableText) joinableI, (JoinableText) joinableJ, parameters)) {
				//joinPointsGroupedIntoJunctions.union(overlap.get(0), overlap.get(1));
					return overlapList;
				}
			} else if ((joinableI instanceof JoinableText && joinableJ.getJoinPoints().size() == 2) || (joinableJ instanceof JoinableText && joinableI.getJoinPoints().size() == 2)) {
				Joinable lineJoinable = (joinableI instanceof JoinableText ? joinableJ : joinableI);
				JoinPoint lineJoinEnd = (overlapList.get(0).getJoinable() instanceof JoinableText ? overlapList.get(1) : overlapList.get(0));
				JoinPoint lineOtherEnd = (lineJoinable.getJoinPoints().get(0) == lineJoinEnd ? lineJoinable.getJoinPoints().get(1) : lineJoinable.getJoinPoints().get(0));
				Line2 line = new Line2(lineOtherEnd.getPoint(), lineJoinEnd.getPoint());
				JoinPoint text = (overlapList.get(0).getJoinable() instanceof JoinableText ? overlapList.get(0) : overlapList.get(1));
				Line2 testLine = new Line2(lineJoinEnd.getPoint(), text.getPoint());
				if (isNumber((SVGText) text.getJoinable().getSVGElement()) || line.isParallelTo(testLine, new Angle(parameters.getTightBondAndTextAngle(), Units.DEGREES))) {
					return overlapList;
				} else {
					text.setRadius(text.getRadius() * parameters.getSmallRadiusExpansion() / parameters.getLargeRadiusExpansion());
					Set<JoinPoint> overlapSet2 = joinableI.overlapWith(joinableJ);
					text.setRadius(text.getRadius() * parameters.getLargeRadiusExpansion() / parameters.getSmallRadiusExpansion());
					if (overlapSet2 != null && line.isParallelTo(testLine, new Angle(parameters.getLooseBondAndTextAngle(), Units.DEGREES))) {
						return overlapList;
					}
				}
			} else {
				return overlapList;
			}
			
			/*if (joinableI instanceof JoinableText && joinableJ instanceof JoinableText) {
				if (doTextsJoin(joinableI, joinableJ)) { 
					joinPointsGroupedIntoJunctions.union(overlap.get(0), overlap.get(1));
				}
			} else {
				joinPointsGroupedIntoJunctions.union(overlap.get(0), overlap.get(1));*/
				/*if (joinableI instanceof HatchedBond) {
					joinables.remove(whichLineIsWhichSingleBond.get(mutuallyExclusivePairs.get(joinableI)));
				}
				if (joinableJ instanceof HatchedBond) {
					joinables.remove(whichLineIsWhichSingleBond.get(mutuallyExclusivePairs.get(joinableJ)));
				}
				if (joinableI instanceof SingleBond) {
					joinables.remove(mutuallyExclusivePairs.inverse().get(joinableI));
				}
				if (joinableJ instanceof SingleBond) {
					joinables.remove(mutuallyExclusivePairs.inverse().get(joinableJ));
				}*/
			//}
			/*if (joinablei instanceof JoinableScriptWord && joinablej instanceof JoinableScriptWord) {
				if (!((JoinableScriptWord) joinablei).getScriptWord().toUnderscoreAndCaretString().equals("H") && !((JoinableScriptWord) joinablej).getScriptWord().toUnderscoreAndCaretString().equals("H")) {
					continue;
				}
			}*/
		}
		return new ArrayList<JoinPoint>();
	}

	private List<JoinPoint> extractAtomLabelsAndGetRemainingJoinPoints(List<Joinable> joinables) {

		List<JoinPoint> remainingJoinPoints = new ArrayList<JoinPoint>();
		Map<Double, List<JoinableText>> listsOfTextsByFontSize = new LinkedHashMap<Double, List<JoinableText>>();
		for (Joinable j : joinables) {
			if (j instanceof JoinableText) {// && isLabel(((JoinableText) j).getSVGElement())) {// && !"1".equals(((JoinableText) j).getSVGElement().getText()) && !"2".equals(((JoinableText) j).getSVGElement().getText()) && !"3".equals(((JoinableText) j).getSVGElement().getText()) && !"4".equals(((JoinableText) j).getSVGElement().getText())) {
				List<JoinableText> joinablesForSize = listsOfTextsByFontSize.get(((JoinableText) j).getSVGElement().getFontSize());
				if (joinablesForSize == null) {
					joinablesForSize = new ArrayList<JoinableText>();
					listsOfTextsByFontSize.put(((JoinableText) j).getSVGElement().getFontSize(), joinablesForSize);
				}
				joinablesForSize.add((JoinableText) j);
			} else {
				remainingJoinPoints.addAll(j.getJoinPoints());
			}
		}
		
		double fontSizeOfLabels = Double.MAX_VALUE;
		list: for (Entry<Double, List<JoinableText>> list : listsOfTextsByFontSize.entrySet()) {
			AreInSameStringDetector sameString = new AreInSameStringDetector(list.getValue(), parameters, false, true);
			ImmutableCollection<Set<Joinable>> groups = sameString.texts.snapshot();
			//List<Integer> labelNumbers = new ArrayList<Integer>();
			Map<Real2Range, Integer> labelNumbers = new LinkedHashMap<Real2Range, Integer>();
			group: for (Set<Joinable> group : groups) {
				List<Joinable> potentialLabelTexts = new ArrayList<Joinable>(group);
				Joinable.sortJoinablesByX(potentialLabelTexts);
				String number = "";
				Real2Range bounds = new Real2Range();
				for (Joinable potentialLabelText : potentialLabelTexts) {
					if (!isLabel((SVGText) potentialLabelText.getSVGElement())) {
						for (Joinable t : group) {
							remainingJoinPoints.addAll(t.getJoinPoints());
						}
						continue group;
					}
					bounds.add(potentialLabelText.getJoinPoints().get(0).getPoint());
					if (isNumber((SVGText) potentialLabelText.getSVGElement())) {
						number += ((SVGText) potentialLabelText.getSVGElement()).getText();
					}
				}
				try {
					labelNumbers.put(bounds, Integer.parseInt(number));
				} catch (NumberFormatException e) {
					
				}
			}
			List<Integer> labelNumbersToBeSorted = new ArrayList<Integer>(labelNumbers.values());
			Collections.sort(labelNumbersToBeSorted);
			int previousPreviousLabel = 0;
			int previousLabel = 0;
			for (Integer i : labelNumbersToBeSorted) {
				if (i - previousLabel > labelNumbersToBeSorted.get(labelNumbersToBeSorted.size() - 1) * parameters.getMaximumLabelSequenceGap()) {
					for (JoinableText t : list.getValue()) {
						remainingJoinPoints.addAll(t.getJoinPoints());
					}
					continue list;
				}
				previousPreviousLabel = previousLabel;
				previousLabel = i;
			}
			if (list.getKey() < fontSizeOfLabels && labelNumbers.size() > 1) {
				if (atomLabelTexts != null) {
					for (JoinableText t : atomLabelTexts) {
						remainingJoinPoints.addAll(t.getJoinPoints());
					}
				}
				atomLabelTexts = list.getValue();
				atomLabelPositionsAndNumbers = labelNumbers;
				fontSizeOfLabels = list.getKey();
			} else {
				for (JoinableText t : list.getValue()) {
					remainingJoinPoints.addAll(t.getJoinPoints());
				}
			}
		}
		
		if (atomLabelTexts == null) {
			atomLabelTexts = new ArrayList<JoinableText>();
			atomLabelPositionsAndNumbers = new HashMap<Real2Range, Integer>();
		}

		return remainingJoinPoints;
	}
	
	private boolean isLetter(SVGText svgElement) {
		return (svgElement.getText() == null ? false : svgElement.getText().matches("[A-Za-z]"));
	}
	
	private boolean isNumber(SVGText svgElement) {
		return (svgElement.getText() == null ? false : svgElement.getText().matches("[0-9]"));
	}
	
	private boolean isLabel(SVGText svgElement) {
		return (svgElement.getText() == null ? false : svgElement.getText().matches("[0-9'a]"));
	}

	protected void createJoinableList() {
		List<Joinable> joinableList = createJoinableList(higherPrimitives.getLineList());
		joinableList.addAll(createJoinableList(derivedPrimitives.getPolygonList()));
		joinableList.addAll(createJoinableList(derivedPrimitives.getPathList()));
		joinableList.addAll(higherPrimitives.getDoubleBondList());
		joinableList.addAll(higherPrimitives.getTripleBondList());
		joinableList.addAll(higherPrimitives.getHatchedBondList());
		joinableList.addAll(higherPrimitives.getLineChargeList());
		//joinableList.addAll(higherPrimitives.getWordList());
		joinableList.addAll(createJoinableList(derivedPrimitives.getTextList()));
		//joinableList.addAll(createJoinableList(derivedPrimitives.getImageList()));
		higherPrimitives.addJoinableList(joinableList);
	}
	
	public List<Joinable> createJoinableList(List<? extends SVGElement> elementList) {
		List<Joinable> joinableList = new ArrayList<Joinable>();
		for (SVGElement element : elementList) {
			List<Joinable> joinables = createJoinables(element);
			joinableList.addAll(joinables);
		}
		return joinableList;
	}

	private List<Joinable> createJoinables(SVGElement element) {
		List<Joinable> joinables = new ArrayList<Joinable>();
		Joinable joinable = null;
		if (element instanceof SVGLine) {
			joinable = new SingleBond(parameters, (SVGLine) element);
			for (AmbiguityHandler.MutuallyExclusiveShortLineTriple triple : ambiguityHandler.mutuallyExclusiveShortLineTriples) {
				if (triple.line == element) {
					triple.singleBond = (SingleBond) joinable;
				}
			}
			for (AmbiguityHandler.MutuallyExclusiveShortLinePairTriple triple : ambiguityHandler.mutuallyExclusiveShortLinePairTriples) {
				if (triple.line1 == element) {
					triple.singleBond1 = (SingleBond) joinable;
				}
				if (triple.line2 == element) {
					triple.singleBond2 = (SingleBond) joinable;
				}
			}
			for (AmbiguityHandler.MutuallyExclusiveShortLineTripleTriple triple : ambiguityHandler.mutuallyExclusiveShortLineTripleTriples) {
				if (triple.line1 == element) {
					triple.singleBond1 = (SingleBond) joinable;
				}
				if (triple.line2 == element) {
					triple.singleBond2 = (SingleBond) joinable;
				}
				if (triple.line3 == element) {
					triple.singleBond3 = (SingleBond) joinable;
				}
			}
			for (AmbiguityHandler.MutuallyExclusiveLinePairPair pair : ambiguityHandler.mutuallyExclusiveLinePairPairs) {
				if (pair.line1 == element) {
					pair.singleBond1 = (SingleBond) joinable;
				}
				if (pair.line2 == element) {
					pair.singleBond2 = (SingleBond) joinable;
				}
			}
			for (AmbiguityHandler.MutuallyExclusiveLineTriplePair pair : ambiguityHandler.mutuallyExclusiveLineTriplePairs) {
				if (pair.line1 == element) {
					pair.singleBond1 = (SingleBond) joinable;
				}
				if (pair.line2 == element) {
					pair.singleBond2 = (SingleBond) joinable;
				}
				if (pair.line3 == element) {
					pair.singleBond3 = (SingleBond) joinable;
				}
			}
		} else if (element instanceof SVGText) {
			if (("+".equals(((SVGText) element).getText()) || "-".equals(((SVGText) element).getText())) && !JoinableText.anyTextsInSameString((SVGText) element, derivedPrimitives.getTextList(), parameters, false, true)) {
				joinable = new Charge(parameters, (SVGText) element);
			} else {
				joinable = new JoinableText(parameters, (SVGText) element);
			}
 		} else if (element instanceof SVGPolygon && ((SVGPolygon) element).createLineList(true).size() == 3) {
 			double shortest = Double.MAX_VALUE;
 			for (SVGLine line : ((SVGPolygon) element).getLineList()) {
 				if (line.getLength() < shortest) {
 					shortest = line.getLength();
 				}
 			}
 			if (shortest > 0.5) {
 				joinable = new WedgeBond(parameters, (SVGPolygon) element);
 			}
 		} else if (element instanceof SVGPolygon && ((SVGPolygon) element).createLineList(true).size() == 4) {
 			Real2Array real2Array = ((SVGPolygon) element).getReal2Array();
			for (int i = 0; i < real2Array.size(); i++) {
 				Real2Array withoutPoint = new Real2Array(real2Array);
 				Real2 deleted = withoutPoint.get(i);
 				withoutPoint.deleteElement(i);
 				SVGPolygon newPoly = new SVGPolygon(withoutPoint);
 				if (newPoly.containsPoint(deleted, 0)) {//withoutPoint.getRange2().includes(((SVGPolygon) element).getReal2Array().get(i))) {
 					if (SimpleBuilder.area(real2Array) / SimpleBuilder.area(newPoly.getReal2Array()) < 0.5) {
 						double shortest = Double.MAX_VALUE;
 						int nearest = -1;
 						for (int p = 0; p < 4; p++) {
 							Real2 point1 = real2Array.get(p);
 							for (Real2 point2 : real2Array) {
 	 							double dist = point1.getDistance(point2);
 	 							if (dist != 0 && dist < shortest) {
 	 								shortest = dist;
 	 								nearest = p;
 	 							}
 	 						}	
 						}
 						SVGPolygon polygon1;
 						SVGPolygon polygon2;
 						if (nearest == 0) {
 							polygon1 = new SVGPolygon(new Real2Array(Arrays.asList(real2Array.get(0), real2Array.get(1), real2Array.get(2))));
 							polygon2 = new SVGPolygon(new Real2Array(Arrays.asList(real2Array.get(0), real2Array.get(3), real2Array.get(2))));
 						} else {
 							polygon1 = new SVGPolygon(new Real2Array(Arrays.asList(real2Array.get(1), real2Array.get(2), real2Array.get(3))));
 							polygon2 = new SVGPolygon(new Real2Array(Arrays.asList(real2Array.get(1), real2Array.get(0), real2Array.get(3))));
 						}	
						WedgeBond wedge1 = new WedgeBond(parameters, polygon1);
						joinables.add(wedge1);
						WedgeBond wedge2 = new WedgeBond(parameters, polygon2);
						joinables.add(wedge2);
 					} else {
	 					((SVGPolygon) element).setReal2Array(withoutPoint);
	 					joinable = new WedgeBond(parameters, (SVGPolygon) element);
	 					break;
 					}
 				}
 			}
 		} else if (element instanceof SVGPath) {
 			try {
 				joinable = new WigglyBond(parameters, (SVGPath) element);
 			} catch (IllegalArgumentException e) {
 				
 			}
 		}
		if (joinable == null && joinables.size() == 0) {
 			LOG.debug("Unknown joinable: " + element);
 		}
		if (joinable != null && joinables.size() == 0) {
			joinables.add(joinable);
		}
		return joinables;
	}

	private void createUnsaturatedBondLists() {
		DoubleBondManager unsaturatedBondManager = new DoubleBondManager(parameters);
		try {
			unsaturatedBondManager.createBondLists(higherPrimitives.getLineList(), startTime + timeout - System.currentTimeMillis());
		} catch (TimeoutException e) {
			throw new UncheckedTimeoutException(e.getMessage());
		}
		//doubleBondManager.removeUsedDoubleBondPrimitives(higherPrimitives.getLineList());
		List<DoubleBond> doubleBondList = unsaturatedBondManager.getDoubleBondList();
		List<TripleBond> tripleBondList = unsaturatedBondManager.getTripleBondList();
		higherPrimitives.setDoubleBondList(doubleBondList);
		higherPrimitives.setTripleBondList(tripleBondList);
		ambiguityHandler.mutuallyExclusiveLinePairPairs = new ArrayList<AmbiguityHandler.MutuallyExclusiveLinePairPair>();
		bond: for (DoubleBond bond : doubleBondList) {
			for (AmbiguityHandler.MutuallyExclusiveShortLinePairTriple pair : ambiguityHandler.mutuallyExclusiveShortLinePairTriples) {
				if ((pair.line1 == bond.getLine(0) && pair.line2 == bond.getLine(1)) || (pair.line1 == bond.getLine(1) && pair.line2 == bond.getLine(0))) {
					pair.doubleBond = bond;
					continue bond;
				}
			}
			ambiguityHandler.mutuallyExclusiveLinePairPairs.add(new AmbiguityHandler.MutuallyExclusiveLinePairPair(bond));
			//higherPrimitives.getLineList().add(bond.getLine(0));
			//higherPrimitives.getLineList().add(bond.getLine(1));
		}
		ambiguityHandler.mutuallyExclusiveLineTriplePairs = new ArrayList<AmbiguityHandler.MutuallyExclusiveLineTriplePair>();
		bond: for (TripleBond bond : tripleBondList) {
			for (AmbiguityHandler.MutuallyExclusiveShortLineTripleTriple triple : ambiguityHandler.mutuallyExclusiveShortLineTripleTriples) {
				if (triple.setIfEquals(bond)) {
					continue bond;
				}
			}
			ambiguityHandler.mutuallyExclusiveLineTriplePairs.add(new AmbiguityHandler.MutuallyExclusiveLineTriplePair(bond));
			//higherPrimitives.getLineList().add(bond.getLine(0));
			//higherPrimitives.getLineList().add(bond.getLine(1));
		}
	}

	/*private void createRawJunctionList() {
		createJoinableList();
		List<Joinable> joinableList = higherPrimitives.getJoinableList();
		List<Junction> rawJunctionList = new ArrayList<Junction>();
		for (int i = 0; i < joinableList.size() - 1; i++) {
			Joinable joinablei = joinableList.get(i);
			for (int j = i + 1; j < joinableList.size(); j++) {
				Joinable joinablej = joinableList.get(j);
				JoinPoint commonPoint = joinablei.getIntersectionPoint(joinablej);
				if (commonPoint != null) {
					Junction junction = new Junction(joinablei, joinablej, commonPoint);
					rawJunctionList.add(junction);
					String junctAttVal = "junct"+"."+rawJunctionList.size();
					junction.addAttribute(new Attribute(SVGElement.ID, junctAttVal));
					if (junction.getCoordinates() == null && commonPoint.getPoint() != null) {
						junction.setCoordinates(commonPoint.getPoint());
					}
					LOG.debug("junct: "+junction.getId()+" between " + joinablei.getClass() + " and " + joinablej.getClass() + " with coords "+junction.getCoordinates()+" "+commonPoint.getPoint());
				}
			}
		}
		higherPrimitives.setRawJunctionList(rawJunctionList);
	}*/

	public SVGElement getSVGRoot() {
		return svgRoot;
	}

	public HigherPrimitives getHigherPrimitives() {
		return higherPrimitives;
	}
	
	public Map<Real2Range, Integer> getAtomLabels() {
		return atomLabelPositionsAndNumbers;
	}

	void draw() {
		draw(new File("target/chem/andy.svg"));
	}
	
	public void draw(File file) {
		SVGG out = drawPrimitivesJoinPointsAndJunctions();
		SVGSVG.wrapAndWriteAsSVG(out, file);
	}

	SVGG drawPrimitivesJoinPointsAndJunctions() {
		SVGG out = new SVGG();
		SVGG circles = new SVGG();
		out.appendChild(circles);
		if (higherPrimitives.getJunctionList() != null) {
			for (Junction j : higherPrimitives.getJunctionList()) {
				Real2 coords = (j.getCoordinates() == null ? new Real2(0, 0) : j.getCoordinates());
				SVGCircle c = new SVGCircle(coords, 1.2);
				c.setFill("#555555");
				c.setOpacity(0.7);
				c.setStrokeWidth(0.0);
				circles.appendChild(c);
				SVGText t = new SVGText(coords.plus(new Real2(1.5, Math.random() * 6)), j.getID());
				circles.appendChild(t);
				for (JoinPoint point : j.getJoinPoints()) {
					SVGLine line = new SVGLine(coords, point.getPoint());
					line.setStrokeWidth(0.05);
					circles.appendChild(line);
				}
			}
		}
		for (SVGText t : getDerivedPrimitives().getTextList()) {
			SVGText o = (SVGText) t.copy();
			out.appendChild(o);
		}
		for (SVGLine l : getDerivedPrimitives().getLineList()) {
			SVGLine o = (SVGLine) l.copy();
			o.setStrokeWidth(0.4);
			out.appendChild(o);
		}
		for (SVGPolygon p : getDerivedPrimitives().getPolygonList()) {
			SVGPolygon o = (SVGPolygon) p.copy();
			o.setStrokeWidth(0.4);
			out.appendChild(o);
		}
		for (SVGPath p : getDerivedPrimitives().getPathList()) {
			SVGPath o = (SVGPath) p.copy();
			o.setStrokeWidth(0.4);
			out.appendChild(o);
		}
		for (Charge t : getHigherPrimitives().getLineChargeList()) {
			if (t.getSVGElement() != null) {
				SVGElement e = (SVGElement) t.getSVGElement().copy();
				out.appendChild(e);
			}
		}
		/*for (SVGImage t : simpleBuilder.getDerivedPrimitives().getImageList()) {
			SVGText e = new SVGText
			out.appendChild(e);
		}*/
		if (getHigherPrimitives().getJoinableList() != null) {
			for (Joinable j : getHigherPrimitives().getJoinableList()) {
				for (JoinPoint p : j.getJoinPoints()) {
					Real2 coords = (p.getPoint() == null ? new Real2(0, 0) : p.getPoint());
					SVGCircle c = new SVGCircle(coords, p.getRadius());
					if (j instanceof SingleBond) {
						c.setFill("#9999FF");
					} else if (j instanceof DoubleBond) {
						c.setFill("#99FF99");
					} else if (j instanceof TripleBond) {
						c.setFill("#CCCCCC");
					} else if (j instanceof HatchedBond) {
						c.setFill("#FF9999");
					} else if (j instanceof WedgeBond) {
						c.setFill("#99FFFF");
					} else if (j instanceof Charge) {
						c.setFill("#FFFF99");
					} else if (j instanceof JoinableText) {
						c.setFill("#FF99FF");
					} else if (j instanceof WigglyBond) {
						c.setFill("#999999");
					}
					c.setOpacity(0.7);
					c.setStrokeWidth(0.0);
					circles.appendChild(c);
					//SVGText t = new SVGText(coords.plus(new Real2(1.5, Math.random() * 6)), j.getId());
					//out.appendChild(t);
				}
			}
		}
		return out;
	}
	
}
