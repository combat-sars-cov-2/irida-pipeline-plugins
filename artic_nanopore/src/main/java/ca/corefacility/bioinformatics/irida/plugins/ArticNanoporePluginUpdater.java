package ca.corefacility.bioinformatics.irida.plugins;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;

import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.google.common.base.Splitter;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ca.corefacility.bioinformatics.irida.exceptions.IridaWorkflowNotFoundException;
import ca.corefacility.bioinformatics.irida.exceptions.PostProcessingException;
import ca.corefacility.bioinformatics.irida.model.sample.Sample;
import ca.corefacility.bioinformatics.irida.model.sample.metadata.MetadataEntry;
import ca.corefacility.bioinformatics.irida.model.sample.metadata.PipelineProvidedMetadataEntry;
import ca.corefacility.bioinformatics.irida.model.workflow.IridaWorkflow;
import ca.corefacility.bioinformatics.irida.model.workflow.analysis.AnalysisOutputFile;
import ca.corefacility.bioinformatics.irida.model.workflow.analysis.type.AnalysisType;
import ca.corefacility.bioinformatics.irida.model.workflow.submission.AnalysisSubmission;
import ca.corefacility.bioinformatics.irida.pipeline.results.updater.AnalysisSampleUpdater;
import ca.corefacility.bioinformatics.irida.service.sample.MetadataTemplateService;
import ca.corefacility.bioinformatics.irida.service.sample.SampleService;
import ca.corefacility.bioinformatics.irida.service.workflow.IridaWorkflowsService;


@JsonIgnoreProperties(ignoreUnknown = true)
class QCReport {
    public float pct_N_bases;
    public float pct_covered_bases;
    public float num_aligned_reads;
	public float longest_no_N_run;
}


class MetadataValue {
    public String header;
    public String value;

    MetadataValue(String header, String value) {
        this.header = header;
        this.value = value;	
    }

}


/**
 * This implements a class used to perform post-processing on the analysis
 * pipeline results to extract information to write into the IRIDA metadata
 * tables. Please see
 * <https://github.com/phac-nml/irida/blob/development/src/main/java/ca/corefacility/bioinformatics/irida/pipeline/results/AnalysisSampleUpdater.java>
 * or the README.md file in this project for more details.
 */
public class ArticNanoporePluginUpdater implements AnalysisSampleUpdater {
	private static final Logger logger = LoggerFactory.getLogger( ArticNanoporePluginUpdater.class);

	private static final String QC_FILE = "qc.json";
	private static final String NEXT_CLADE_FILE = "nextclade.tsv";
	private static final String PANGOLIN_FILE = "pangolin.csv";

	private static final Splitter SPLITTER = Splitter.on('\t');
	
	private final MetadataTemplateService metadataTemplateService;
	private final SampleService sampleService;
	private final IridaWorkflowsService iridaWorkflowsService;


	/**
	 * Builds a new {@link ArticNanoporePluginUpdater} with the given services.
	 * 
	 * @param metadataTemplateService The metadata template service.
	 * @param sampleService           The sample service.
	 * @param iridaWorkflowsService   The irida workflows service.
	 */
	public ArticNanoporePluginUpdater(MetadataTemplateService metadataTemplateService, SampleService sampleService,
			IridaWorkflowsService iridaWorkflowsService) {
		this.metadataTemplateService = metadataTemplateService;
		this.sampleService = sampleService;
		this.iridaWorkflowsService = iridaWorkflowsService;
	}
	
	/**
	 * Parses a line of the results file and gets a Map linking the column to the
	 * value in the line. (e.g., "N50 value" => "100").
	 * 
	 * @param columnNames        A List of names of the columns in the results file.
	 * @param line               The line to parse.
	 * @param singleColumnPrefix A prefix for a special case in the staramr results
	 *                           where the column prefix is constant but the suffix
	 *                           changes. Set to null to ignore.
	 * @param resultsFile        The specific file being parsed (for error
	 *                           messages).
	 * @param analysis           The analysis submission being parsed (for error
	 *                           messages).
	 * @return A Map linking the column to the value for the line.
	 * @throws PostProcessingException If there was an error parsing the results.
	 */
	
	private Map<String, String> getDataMapForLine(List<String> columnNames, String line,
			Path resultsFile, AnalysisSubmission analysis) throws PostProcessingException {
		Map<String, String> dataMap = new HashMap<>();

		List<String> values = SPLITTER.splitToList(line);

		if (columnNames.size() != values.size()) {
			throw new PostProcessingException("Mismatch in number of column names [" + columnNames.size()
					+ "] and number of files [" + values.size() + "] in results file [" + resultsFile + "]");
		}

		for (int i = 0; i < columnNames.size(); i++) {
			dataMap.put(columnNames.get(i), values.get(i));
		}
		return dataMap;
	}
	
	/**
	 * Gets the qc results from the given output file.
	 * 
	 * @param qcFilePath The QC output file containing the results.
	 * @param analysis        The {@link AnalysisSubmission} containing the results.
	 * @return A {@link Map} storing the results from staramr, keyed by the metadata
	 *         field name.
	 * @throws IOException             If there was an issue reading the file.
	 * @throws PostProcessingException If there was an issue parsing the file.
	 */
	private Map<String, PipelineProvidedMetadataEntry> getQCResults(Path qcFilePath,
			AnalysisSubmission analysis) throws IOException, PostProcessingException {
		
		ObjectMapper mapper = new ObjectMapper();
		Map<String, PipelineProvidedMetadataEntry> results = new HashMap<>();
		
		@SuppressWarnings("resource")
		String scannedFile = new Scanner(new BufferedReader(new FileReader(qcFilePath.toFile()))).useDelimiter("\\Z").next();
		QCReport qc_results = mapper.readValue(scannedFile, QCReport.class);

		results.put("pct_N_bases", new PipelineProvidedMetadataEntry(Float.toString(qc_results.pct_N_bases), "Percentage N Bases", analysis));
		results.put("pct_covered_bases", new PipelineProvidedMetadataEntry(Float.toString(qc_results.pct_covered_bases), "Percentage Covered Bases", analysis));
		results.put("num_aligned_reads", new PipelineProvidedMetadataEntry(Float.toString(qc_results.num_aligned_reads), "Number of Aligned Reads", analysis));
		results.put("longest_no_N_run", new PipelineProvidedMetadataEntry(Float.toString(qc_results.longest_no_N_run), "Longest of Runs", analysis));
		
		return results;
	}

	/**
	 * Gets the nextclade results from the given output file.
	 * 
	 * @param nextcladeFilePath The nextclade output file containing the results.
	 * @param analysis        The {@link AnalysisSubmission} containing the results.
	 * @return A {@link Map} storing the results from staramr, keyed by the metadata
	 *         field name.
	 * @throws IOException             If there was an issue reading the file.
	 * @throws PostProcessingException If there was an issue parsing the file.
	 */
	private Map<String, PipelineProvidedMetadataEntry> getNextCladeResults(Path nextcladeFilePath,
			AnalysisSubmission analysis) throws IOException, PostProcessingException {
		final int MIN_TOKENS = 2;

		Map<String, PipelineProvidedMetadataEntry> results = new HashMap<>();
		Map<String, String> dataMap;

		@SuppressWarnings("resource")
		BufferedReader reader = new BufferedReader(new FileReader(nextcladeFilePath.toFile()));
		String line = reader.readLine();
		List<String> columnNames = SPLITTER.splitToList(line);
		
		if (columnNames.size() < MIN_TOKENS) {
			throw new PostProcessingException("Invalid number of columns in nextclade results file [" + nextcladeFilePath
					+ "], expected at least [" + MIN_TOKENS + "] got [" + columnNames.size() + "]");
		}

		line = reader.readLine();
		if (line == null || line.length() == 0) {
			dataMap = new HashMap<>();
			logger.info(
					"Got empty results for nextclade file [" + nextcladeFilePath + "] for analysis submission " + analysis);
		} else {
			dataMap = getDataMapForLine(columnNames, line, nextcladeFilePath, analysis);
		}
		results.put("aaSubstitutions", new PipelineProvidedMetadataEntry(dataMap.get("aaSubstitutions"), "Amino Acids Substitutions", analysis));
		results.put("substitutions", new PipelineProvidedMetadataEntry(dataMap.get("substitutions"), "Variants", analysis));
		results.put("clade", new PipelineProvidedMetadataEntry(dataMap.get("clade"), "Clades", analysis));
		
		line = reader.readLine();
		
		if (line == null) {
			return results;
		} else {
			throw new PostProcessingException("Invalid number of results in nextclade results file [" + nextcladeFilePath
					+ "], expected only one line of results but got multiple lines");
		}
	}

/**
	 * Gets the pangolin results from the given output file.
	 * 
	 * @param pangolinFilePath The nextclade output file containing the results.
	 * @param analysis        The {@link AnalysisSubmission} containing the results.
	 * @return A {@link Map} storing the results from staramr, keyed by the metadata
	 *         field name.
	 * @throws IOException             If there was an issue reading the file.
	 * @throws PostProcessingException If there was an issue parsing the file.
	 */
	private Map<String, PipelineProvidedMetadataEntry> getPangolinResults(Path pangolinFilePath,
			AnalysisSubmission analysis) throws IOException, PostProcessingException {
		final int MIN_TOKENS = 2;

		Map<String, PipelineProvidedMetadataEntry> results = new HashMap<>();
		Map<String, String> dataMap;

		@SuppressWarnings("resource")
		BufferedReader reader = new BufferedReader(new FileReader(pangolinFilePath.toFile()));
		String line = reader.readLine();
		List<String> columnNames = SPLITTER.splitToList(line);
		
		if (columnNames.size() < MIN_TOKENS) {
			throw new PostProcessingException("Invalid number of columns in pangolin results file [" + pangolinFilePath
					+ "], expected at least [" + MIN_TOKENS + "] got [" + columnNames.size() + "]");
		}

		line = reader.readLine();
		if (line == null || line.length() == 0) {
			dataMap = new HashMap<>();
			logger.info(
					"Got empty results for pangolin file [" + pangolinFilePath + "] for analysis submission " + analysis);
		} else {
			dataMap = getDataMapForLine(columnNames, line, pangolinFilePath, analysis);
		}
		logger.debug("# DataMap: " + dataMap);
		results.put("lineage", new PipelineProvidedMetadataEntry(dataMap.get("lineage"), "Lineage", analysis));
		logger.debug("# Results: " + results);

		line = reader.readLine();
		
		if (line == null) {
			return results;
		} else {
			throw new PostProcessingException("Invalid number of results in nextclade results file [" + pangolinFilePath
					+ "], expected only one line of results but got multiple lines");
		}
	}


	/**
	 * Code to perform the actual update of the {@link Sample}s passed in the
	 * collection.
	 * 
	 * @param samples  A collection of {@link Sample}s that were passed to this
	 *                 pipeline.
	 * @param analysis The {@link AnalysisSubmission} object corresponding to this
	 *                 analysis pipeline.
	 */
	@Override
	public void update(Collection<Sample> samples, AnalysisSubmission analysis) throws PostProcessingException {
		if (samples == null) {
			throw new IllegalArgumentException("samples is null");
		} else if (analysis == null) {
			throw new IllegalArgumentException("analysis is null");
		} else if (samples.size() != 1) {
			// In this particular pipeline, only one sample should be run at a time so I
			// verify that the collection of samples I get has only 1 sample
			throw new IllegalArgumentException(
					"samples size=" + samples.size() + " is not 1 for analysisSubmission=" + analysis.getId());
		}

		// extract the 1 and only sample (if more than 1, would have thrown an exception
		// above)
		final Sample sample = samples.iterator().next();

		// extracts paths to the analysis result files
		AnalysisOutputFile qcFile = analysis.getAnalysis().getAnalysisOutputFile(QC_FILE);
		AnalysisOutputFile nextCladeFile = analysis.getAnalysis().getAnalysisOutputFile(NEXT_CLADE_FILE);
		AnalysisOutputFile pangolinFile = analysis.getAnalysis().getAnalysisOutputFile(PANGOLIN_FILE);
		
		Path qcFilePath = qcFile.getFile();
		Path nextCladeFilePath = nextCladeFile.getFile();
		Path pangolinFilePath = pangolinFile.getFile(); 
		
		Map<String, MetadataEntry> stringEntries = new HashMap<>();
		try {
			IridaWorkflow iridaWorkflow = iridaWorkflowsService.getIridaWorkflow(analysis.getWorkflowId());
			String workflowVersion = iridaWorkflow.getWorkflowDescription().getVersion();

			Map<String, PipelineProvidedMetadataEntry> qcResult = getQCResults(qcFilePath, analysis);
			Map<String, PipelineProvidedMetadataEntry> nextcladeResult = getNextCladeResults(nextCladeFilePath, analysis);
			Map<String, PipelineProvidedMetadataEntry> pangolinResult = getPangolinResults(pangolinFilePath, analysis);

			List<Map<String, PipelineProvidedMetadataEntry>> resultsList = Arrays.asList(qcResult, nextcladeResult, pangolinResult);
			logger.debug("# Result List" + resultsList);

			for (Map<String, PipelineProvidedMetadataEntry> result : resultsList) {
				for (String fieldName : result.keySet()) {
					stringEntries.put(appendVersion(result.get(fieldName).getType(), workflowVersion), result.get(fieldName));
				}
			}

			Set<MetadataEntry> metadataSet = metadataTemplateService.convertMetadataStringsToSet(stringEntries);
			sampleService.mergeSampleMetadata(sample, metadataSet);
			
		} catch (IOException e) {
			logger.error("Got IOException", e);
			throw new PostProcessingException("Error parsing JSON from results", e);		
		} catch (IridaWorkflowNotFoundException e) {
			throw new PostProcessingException("Workflow not found, for id=" + analysis.getWorkflowId(), e);
		} catch (Exception e) {
			logger.error("Got Exception", e);
			throw e;
		}
	}
	/**
	 * Appends the name and version together for a metadata field name.
	 * 
	 * @param name    The name.
	 * @param version The version.
	 * @return The appended name and version.
	 */
	private String appendVersion(String name, String version) {
		return name + " / " + version;
	}

	/**
	 * The {@link AnalysisType} this {@link AnalysisSampleUpdater} corresponds to.
	 * 
	 * @return The {@link AnalysisType} this {@link AnalysisSampleUpdater}
	 *         corresponds to.
	 */
	@Override
	public AnalysisType getAnalysisType() {
		return ArticNanoporePlugin.DEFAULT;
	}
}
