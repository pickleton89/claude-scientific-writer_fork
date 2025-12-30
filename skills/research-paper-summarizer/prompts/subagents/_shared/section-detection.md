# Section Detection Utility

You are analyzing a scientific paper to detect section boundaries. Your job is to identify the page numbers where each major section begins and ends.

## Input

<pdf_content>
{{PDF_FIRST_PAGES}}
</pdf_content>

<total_pages>{{TOTAL_PAGES}}</total_pages>

## Your Task

Scan the provided PDF content and identify the page numbers where each section begins. Look for:

1. **Section headers** - Bold/capitalized text like "ABSTRACT", "INTRODUCTION", "METHODS", "RESULTS", "DISCUSSION"
2. **Common variants** - "Materials and Methods", "Experimental Procedures", "Methodology", "Findings", "Conclusions"
3. **Visual cues** - Section numbers, horizontal rules, significant whitespace before headers

## Section Keywords to Detect

| Section | Keywords to Look For |
|---------|---------------------|
| Abstract | "Abstract", "Summary" |
| Introduction | "Introduction", "Background" |
| Methods | "Methods", "Materials and Methods", "Experimental Procedures", "Experimental Section", "Methodology" |
| Results | "Results", "Findings" |
| Discussion | "Discussion", "Conclusions", "Concluding Remarks" |
| References | "References", "Bibliography", "Literature Cited" |

## Output Format

Return a JSON object with the detected section boundaries:

```json
{
  "detected_sections": {
    "abstract": {"start": 1, "end": 1, "confidence": "high"},
    "introduction": {"start": 1, "end": 3, "confidence": "high"},
    "methods": {"start": 4, "end": 7, "confidence": "medium"},
    "results": {"start": 8, "end": 14, "confidence": "medium"},
    "discussion": {"start": 15, "end": 18, "confidence": "high"},
    "references": {"start": 19, "end": 22, "confidence": "high"}
  },
  "detection_notes": "Methods section header found on page 4. Results not explicitly labeled but inferred from content transition.",
  "fallback_recommended": false
}
```

## Confidence Levels

- **high**: Explicit section header found
- **medium**: Section inferred from content/structure
- **low**: Estimated based on page ratios

## Fallback Estimation Rules

If explicit section headers are not found, use these page ratio estimates:

| Section | Page Range (% of total) |
|---------|------------------------|
| Abstract + Introduction | Pages 1-3 or first 15% |
| Methods | 15-30% of total pages |
| Results | 30-65% of total pages |
| Discussion | 65-85% of total pages |
| References | 85-100% of total pages |

Set `"fallback_recommended": true` if you had to rely heavily on estimation.

## Important Notes

- Some papers combine Results and Discussion into one section
- Review articles may not have a distinct Methods section
- Supplementary materials are typically after References - do not include in section mapping
- If a section spans less than 1 page, use the same start and end page number

Begin your analysis now. Output ONLY the JSON object, no additional text.
