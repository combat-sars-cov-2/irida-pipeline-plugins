package ca.corefacility.bioinformatics.irida.plugins;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Collection;
import java.util.HashMap;
import static java.util.Map.entry;
import java.util.Map;
import java.util.Set;
import java.util.Scanner;

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


import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.MappingIterator;
import com.fasterxml.jackson.databind.ObjectReader;
import com.fasterxml.jackson.dataformat.csv.CsvMapper;
import com.fasterxml.jackson.dataformat.csv.CsvSchema;

@JsonIgnoreProperties(ignoreUnknown = true)
class QCReport {
    public float pct_N_bases;
    public float pct_covered_bases;
    public float num_aligned_reads;
	public float longest_no_N_run;
}

class NextClade {
	public String aaSubstitutions;
	public String substitutions;
	public String clade;
}

class Pangolin {
	public String lineage; //B.1.177
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
public class ArticIlluminaPluginUpdater implements AnalysisSampleUpdater {
	//private static final Logger logger = LoggerFactory.getLogger(ArticIlluminaPluginUpdater.class);

	private static final String QC_FILE = "qc.json";
	private static final String NEXT_CLADE_FILE = "nextclade.tsv";
	private static final String PANGOLIN_FILE = "pangolin.csv";


	private final MetadataTemplateService metadataTemplateService;
	private final SampleService sampleService;
	private final IridaWorkflowsService iridaWorkflowsService;

	/**
	 * Builds a new {@link ArticIlluminaPluginUpdater} with the given services.
	 * 
	 * @param metadataTemplateService The metadata template service.
	 * @param sampleService           The sample service.
	 * @param iridaWorkflowsService   The irida workflows service.
	 */
	public ArticIlluminaPluginUpdater(MetadataTemplateService metadataTemplateService, SampleService sampleService,
			IridaWorkflowsService iridaWorkflowsService) {
		this.metadataTemplateService = metadataTemplateService;
		this.sampleService = sampleService;
		this.iridaWorkflowsService = iridaWorkflowsService;
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

		try {
			Map<String, MetadataEntry> metadataEntries = new HashMap<>();

			// get information about the workflow (e.g., version and name)
			IridaWorkflow iridaWorkflow = iridaWorkflowsService.getIridaWorkflow(analysis.getWorkflowId());
			String workflowVersion = iridaWorkflow.getWorkflowDescription().getVersion();
			String workflowName = iridaWorkflow.getWorkflowDescription().getName();

			// Metadata: ETL file list (Hashmap)
			HashMap<Path, String[]> etlFiles = new HashMap<Path, String[]>();
			etlFiles.put(qcFilePath, new String[] {"pct_N_bases", "pct_covered_bases", "num_aligned_reads", "longest_no_N_run"});
			etlFiles.put(nextCladeFilePath, new String [] {"aa_substitutions", "variants", "clades"});
			etlFiles.put(pangolinFilePath, new String [] {"lineage"});

			CsvMapper m = new CsvMapper();
			// load the metadata
			for (Path i : etlFiles.keySet()) {
				@SuppressWarnings("resource")
				// map the results into a Map
				ObjectMapper mapper = new ObjectMapper();
				Map<String, MetadataValue> metaFields = new HashMap<>();
				switch(i.getFileName().toString())
				{
					case "qc.json":
						// QC:  Read the outputfile
						String scannedFile = new Scanner(new BufferedReader(new FileReader(i.toFile()))).useDelimiter("\\Z").next();
						QCReport qc_results = mapper.readValue(scannedFile, QCReport.class);
						metaFields = Map.ofEntries(
								entry("pct_N_bases", new MetadataValue("Percentage N Bases", Float.toString(qc_results.pct_N_bases))),
								entry("pct_covered_bases", new MetadataValue("Percentage Covered Bases", Float.toString(qc_results.pct_covered_bases))),
								entry("num_aligned_reads", new MetadataValue("Number of Aligned Reads", Float.toString(qc_results.num_aligned_reads))),
								entry("longest_no_N_run", new MetadataValue("Longest of Runs", Float.toString(qc_results.longest_no_N_run))));
						break;
					case "nextclade.tsv":
						CsvSchema n_schema = m.schemaFor(NextClade.class).withColumnSeparator('\t');
						MappingIterator<Map<String, String>> tsv_results = mapper
								.readerForMapOf(NextClade.class)
								.with(n_schema)
								.readValues(i.toFile());
						Map<String, String> tsv_row = tsv_results.nextValue();
						metaFields = Map.ofEntries(
								entry("aaSubstitutions", new MetadataValue("Amino Acids Substitutions", tsv_row.get("aaSubstitutions"))),
								entry("substitutions", new MetadataValue("Variants", tsv_row.get("substitutions"))),
								entry("clades", new MetadataValue("Clades", tsv_row.get("clades"))));
						break;
					case "pangolin.csv":
						CsvSchema p_schema = m.schemaFor(Pangolin.class).withColumnSeparator('\t');
						/*CsvSchema schema = CsvSchema.builder()
								.addColumn("lineage")
								.build();*/
						MappingIterator<Map<String, String>> csv_results = mapper
								.readerForMapOf(Pangolin.class)
								.with(p_schema)
								.readValues(i.toFile());
						Map<String, String> csv_row = csv_results.nextValue();
						metaFields = Map.ofEntries(
								entry("lineage", new MetadataValue("Lineage", csv_row.get("lineage"))));
						break;
				}
				// @formatter:on */
				metaFields.put("workflow_name", new MetadataValue("Workflow Name", workflowName));
				metaFields.put("workflow_version", new MetadataValue("Workflow Version", workflowVersion));

				metaFields.entrySet().forEach(entry -> {
					PipelineProvidedMetadataEntry metadataEntry =
							new PipelineProvidedMetadataEntry(entry.getValue().value, "text", analysis);
					metadataEntries.put(entry.getValue().header, metadataEntry);
				});

				Set<MetadataEntry> metadataSet = metadataTemplateService.convertMetadataStringsToSet(metadataEntries);
				samples.forEach(s -> {
					sampleService.mergeSampleMetadata(s, metadataSet);
				});
			}
		} catch (JsonProcessingException e) {
			throw new PostProcessingException("Error parsing JSON from results", e);
		} catch (IOException e) {
			throw new PostProcessingException("Error parsing hash file", e);
		} catch (IridaWorkflowNotFoundException e) {
			throw new PostProcessingException("Could not find workflow for id=" + analysis.getWorkflowId(), e);
		}
	}

	/**
	 * The {@link AnalysisType} this {@link AnalysisSampleUpdater} corresponds to.
	 * 
	 * @return The {@link AnalysisType} this {@link AnalysisSampleUpdater}
	 *         corresponds to.
	 */
	@Override
	public AnalysisType getAnalysisType() {
		return ArticIlluminaPlugin.DEFAULT;
	}
}
