# Instalar pacotes se não tiver
#if(!requireNamespace("Biostrings", quietly = TRUE)) install.packages("BiocManager")
#if(!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
#if(!requireNamespace("ggforce", quietly = TRUE)) install.packages("ggforce")
#if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("Biostrings")\

library(Biostrings)
library(ggplot2)
library(ggforce)

# -----------------------------
# 1. Ler arquivos FASTA
# -----------------------------
ref_fasta <- "BG1_03444.fasta"       # caminho do seu arquivo de referência
primers_fasta <- "primers.fasta" # arquivo contendo primers (forward e reverse)

# Ler sequência de referência
ref_seq <- readDNAStringSet(ref_fasta)
ref_seq <- ref_seq[[1]]  # pegar a primeira sequência se houver mais de uma

# Ler primers
primers <- readDNAStringSet(primers_fasta)

# -----------------------------
# 2. Função para localizar primers
# -----------------------------
find_primers <- function(ref, primer_set){
  results <- data.frame(
    Primer = character(),
    Start = numeric(),
    End = numeric(),
    Strand = character(),
    Score = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(i in seq_along(primer_set)){
    primer <- primer_set[[i]]
    
    # Forward
    fwd <- pairwiseAlignment(primer, ref, type = "local")
    results <- rbind(results, data.frame(
      Primer = names(primer_set)[i],
      Start = start(subject(fwd)),
      End = end(subject(fwd)),
      Strand = "+",
      Score = score(fwd)
    ))
    
    # Reverse (buscar no reverse complement)
    rev_primer <- reverseComplement(primer)
    rev <- pairwiseAlignment(rev_primer, ref, type = "local")
    results <- rbind(results, data.frame(
      Primer = names(primer_set)[i],
      Start = start(subject(rev)),
      End = end(subject(rev)),
      Strand = "-",
      Score = score(rev)
    ))
  }
  return(results)
}

# -----------------------------
# 3. Rodar validação
# -----------------------------
validation <- find_primers(ref_seq, primers)
print(validation)

# -----------------------------
# 4. Plotar localização dos primers
# -----------------------------
library(ggplot2)
library(ggforce)

# Criar coluna para posicionamento visual das barras (direção)
validation$ymin <- as.numeric(factor(validation$Primer)) - 0.2
validation$ymax <- as.numeric(factor(validation$Primer)) + 0.2

# Adicionar setas para indicar direção
validation$arrow_dir <- ifelse(validation$Strand == "+", 1, -1)

ggplot() +
  # Linha base da sequência
  geom_segment(aes(x = 1, xend = max(validation$End)+5, y = 0, yend = 0), color="black", size=0.5) +
  
  # Barras dos primers
  geom_rect(data=validation, aes(xmin=Start, xmax=End, ymin=ymin, ymax=ymax, fill=Strand), color="black") +
  
  # Setas indicando direção
  geom_segment(data=validation, aes(x=Start, xend=End, y=(ymin+ymax)/2, yend=(ymin+ymax)/2),
               arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "closed"),
               color="black") +
  
  # Score dentro da barra
  geom_text(data=validation, aes(x=(Start+End)/2, y=(ymin+ymax)/2, label=paste0("Score=", round(Score,1))),
            color="white", size=3) +
  
  # Escala de cores
  scale_fill_manual(values = c("+" = "#1f78b4", "-" = "#e31a1c")) +
  
  # Labels
  labs(title = "Localização dos Primers na Sequência de Referência (BG1_03444)",
       x = "Posição na referência (nt)",
       y = "Primer",
       fill = "Sentido") +
  
  # Ajustes de tema
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=10))
