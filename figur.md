graph TD
    %% Girdiler
    A[(VCF File)] --> C(Data Parsing & SNP Filtering)
    B[(Reference Genome .fa.gz)] --> C
    
    %% İşlem Adımları
    C -->|Biallelic SNPs| D(Flanking Sequence Extraction)
    D --> E(Allele Sequence Generation: REF vs ALT)
    
    E --> F{Restriction Enzyme Screening}
    
    %% Karar Ağacı
    F -->|No Polymorphic Cut| G[Discard / No CAPS Marker]
    F -->|Polymorphic Cut Found| H(Primer Design via Primer3)
    
    H -->|Optimization: Tm, GC, Size| I[Successful CAPS/TaqMan Marker]
    
    %% Çıktılar
    I --> J{Results & Visualization}
    
    J --> K[In Silico Gel Simulator]
    J --> L[Physical Chromosome Map]
    J --> M[\Excel Report Export/]
    
    %% Stil Ayarları (Siberpunk/Koyu Tema Uyumlu)
    style A fill:#1e293b,stroke:#0ea5e9,stroke-width:2px,color:#fff
    style B fill:#1e293b,stroke:#0ea5e9,stroke-width:2px,color:#fff
    style C fill:#0f172a,stroke:#10b981,stroke-width:2px,color:#fff
    style D fill:#0f172a,stroke:#10b981,stroke-width:2px,color:#fff
    style E fill:#0f172a,stroke:#10b981,stroke-width:2px,color:#fff
    style F fill:#334155,stroke:#f59e0b,stroke-width:2px,color:#fff
    style G fill:#7f1d1d,stroke:#ef4444,stroke-width:2px,color:#fff
    style H fill:#0f172a,stroke:#10b981,stroke-width:2px,color:#fff
    style I fill:#064e3b,stroke:#10b981,stroke-width:3px,color:#fff
    style J fill:#1e293b,stroke:#8b5cf6,stroke-width:2px,color:#fff
    style K fill:#312e81,stroke:#8b5cf6,stroke-width:2px,color:#fff
    style L fill:#312e81,stroke:#8b5cf6,stroke-width:2px,color:#fff
    style M fill:#064e3b,stroke:#10b981,stroke-width:2px,color:#fff