## Installation 

Install the fastRUVIII R package with devtools: 

```r
devtools::install_github('Neal-Liu1/fastRUVIII')
```

```mermaid
flowchart TD
  A[Raw counts<br/>(multi-modal)]

  A --> B{Known biology?}
  B -->|Yes| P[Create PRPC()]
  B -->|No| C[Clustering<br/>(after basic normalization)]
  C --> N[Find Corrected /<br/> Multimodal Neighbours]
  N --> P

  A --> R[Find NCGs]
  R --> S[Assess performance]
  S --> F[fastRUVIII]

  P --> F
  F --> W[Assess W]
  F --> D[Downstream analysis]
```
