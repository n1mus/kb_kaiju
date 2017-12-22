
package us.kbase.kbkaiju;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: KaijuInputParams</p>
 * <pre>
 * Kaiju App Input Params
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "reads_ref",
    "tax_levels",
    "db_type",
    "seg_filter",
    "greedy_run_mode",
    "min_match_length",
    "greedy_min_match_score",
    "greedy_allowed_mismatches"
})
public class KaijuInputParams {

    @JsonProperty("workspace_name")
    private java.lang.String workspaceName;
    @JsonProperty("reads_ref")
    private java.lang.String readsRef;
    @JsonProperty("tax_levels")
    private List<String> taxLevels;
    @JsonProperty("db_type")
    private java.lang.String dbType;
    @JsonProperty("seg_filter")
    private Long segFilter;
    @JsonProperty("greedy_run_mode")
    private Long greedyRunMode;
    @JsonProperty("min_match_length")
    private Long minMatchLength;
    @JsonProperty("greedy_min_match_score")
    private Double greedyMinMatchScore;
    @JsonProperty("greedy_allowed_mismatches")
    private Long greedyAllowedMismatches;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("workspace_name")
    public java.lang.String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(java.lang.String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public KaijuInputParams withWorkspaceName(java.lang.String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("reads_ref")
    public java.lang.String getReadsRef() {
        return readsRef;
    }

    @JsonProperty("reads_ref")
    public void setReadsRef(java.lang.String readsRef) {
        this.readsRef = readsRef;
    }

    public KaijuInputParams withReadsRef(java.lang.String readsRef) {
        this.readsRef = readsRef;
        return this;
    }

    @JsonProperty("tax_levels")
    public List<String> getTaxLevels() {
        return taxLevels;
    }

    @JsonProperty("tax_levels")
    public void setTaxLevels(List<String> taxLevels) {
        this.taxLevels = taxLevels;
    }

    public KaijuInputParams withTaxLevels(List<String> taxLevels) {
        this.taxLevels = taxLevels;
        return this;
    }

    @JsonProperty("db_type")
    public java.lang.String getDbType() {
        return dbType;
    }

    @JsonProperty("db_type")
    public void setDbType(java.lang.String dbType) {
        this.dbType = dbType;
    }

    public KaijuInputParams withDbType(java.lang.String dbType) {
        this.dbType = dbType;
        return this;
    }

    @JsonProperty("seg_filter")
    public Long getSegFilter() {
        return segFilter;
    }

    @JsonProperty("seg_filter")
    public void setSegFilter(Long segFilter) {
        this.segFilter = segFilter;
    }

    public KaijuInputParams withSegFilter(Long segFilter) {
        this.segFilter = segFilter;
        return this;
    }

    @JsonProperty("greedy_run_mode")
    public Long getGreedyRunMode() {
        return greedyRunMode;
    }

    @JsonProperty("greedy_run_mode")
    public void setGreedyRunMode(Long greedyRunMode) {
        this.greedyRunMode = greedyRunMode;
    }

    public KaijuInputParams withGreedyRunMode(Long greedyRunMode) {
        this.greedyRunMode = greedyRunMode;
        return this;
    }

    @JsonProperty("min_match_length")
    public Long getMinMatchLength() {
        return minMatchLength;
    }

    @JsonProperty("min_match_length")
    public void setMinMatchLength(Long minMatchLength) {
        this.minMatchLength = minMatchLength;
    }

    public KaijuInputParams withMinMatchLength(Long minMatchLength) {
        this.minMatchLength = minMatchLength;
        return this;
    }

    @JsonProperty("greedy_min_match_score")
    public Double getGreedyMinMatchScore() {
        return greedyMinMatchScore;
    }

    @JsonProperty("greedy_min_match_score")
    public void setGreedyMinMatchScore(Double greedyMinMatchScore) {
        this.greedyMinMatchScore = greedyMinMatchScore;
    }

    public KaijuInputParams withGreedyMinMatchScore(Double greedyMinMatchScore) {
        this.greedyMinMatchScore = greedyMinMatchScore;
        return this;
    }

    @JsonProperty("greedy_allowed_mismatches")
    public Long getGreedyAllowedMismatches() {
        return greedyAllowedMismatches;
    }

    @JsonProperty("greedy_allowed_mismatches")
    public void setGreedyAllowedMismatches(Long greedyAllowedMismatches) {
        this.greedyAllowedMismatches = greedyAllowedMismatches;
    }

    public KaijuInputParams withGreedyAllowedMismatches(Long greedyAllowedMismatches) {
        this.greedyAllowedMismatches = greedyAllowedMismatches;
        return this;
    }

    @JsonAnyGetter
    public Map<java.lang.String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(java.lang.String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public java.lang.String toString() {
        return ((((((((((((((((((((("KaijuInputParams"+" [workspaceName=")+ workspaceName)+", readsRef=")+ readsRef)+", taxLevels=")+ taxLevels)+", dbType=")+ dbType)+", segFilter=")+ segFilter)+", greedyRunMode=")+ greedyRunMode)+", minMatchLength=")+ minMatchLength)+", greedyMinMatchScore=")+ greedyMinMatchScore)+", greedyAllowedMismatches=")+ greedyAllowedMismatches)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
