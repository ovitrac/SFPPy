#!/usr/bin/env python3
# example5.py
#
# ---------------------------------------------
# Knowledge-Based Question Answering using RAG
# Integrated into SFPPy / Generative Simulation Initiative ✨
# Maintainer: olivier.vitrac@gmail.com
# Created: 2025-07-25
# License: MIT
# ---------------------------------------------

"""
§ Example 5: Retrieval-Augmented Generation (RAG) over Legal Knowledge Base (Markdown)

Usage:
    python example5.py                      # loads existing index and prompts user
    python example5.py --rebuild           # rebuilds the index then prompts user
    python example5.py --rebuild "query"    # rebuilds and runs specific query
    python example5.py "query"              # loads and runs specific query

This script demonstrates how to perform a RAG (Retrieval-Augmented Generation) query on a local Markdown Knowledge Base using the Mistral model served through Ollama, and LlamaIndex for indexing and querying.

---
§ Why Mistral and Markdown KB?
- **Mistral** is fast, instruction-tuned, open-source, and fits on 8 GB VRAM.
- **Markdown** is a natural format for scientific and legal documentation (e.g. EU Regulation 10/2011).
- The KB is re-indexed only if files change, enabling prefetch for performance and auditability.

---
§ Auditability:
- Vector index is persistent, reproducible, and versioned
- KB source files are preserved under `docs/KB/`
- Query and response are printed with traceable file origins

---
§ Dependencies (non-SFPPy):
To install, activate your environment (e.g. `conda activate sfppy`), then:

    pip install llama-index llama-index-embeddings-huggingface \
                llama-index-llms-ollama llama-index-vector-stores-chroma \
                llama-index-vector-stores sentence-transformers chromadb

Optional:
    pip install unstructured  # for advanced HTML/PDF ingestion (not used here)

🌱 Part of the Generative Simulation Initiative – https://github.com/ovitrac/SFPPy
Contact: olivier.vitrac@gmail.com
"""

# ---------------------------------------------
#                CONFIGURATION
# ---------------------------------------------

import os
import sys
import time
import argparse
from pathlib import Path
from llama_index.core import SimpleDirectoryReader, VectorStoreIndex, StorageContext, load_index_from_storage
from llama_index.embeddings.huggingface import HuggingFaceEmbedding
from llama_index.llms.ollama import Ollama
from llama_index.vector_stores.chroma import ChromaVectorStore
import chromadb

# Default question
question_default = "Question: What is the definition of a functional barrier in EU Regulation 10/2011?"
question_default = "What does it mean migration?"

# Paths
mainfolder = Path(os.getcwd())
kbdir = mainfolder / "docs" / "KB"
indexdir = mainfolder / "docs" / "KBIndex"
indexdir.mkdir(parents=True, exist_ok=True)

# CLI args
parser = argparse.ArgumentParser(description="SFPPy RAG over Markdown KB")
parser.add_argument("query", nargs="?", help="Question to ask (optional)")
parser.add_argument("--rebuild", action="store_true", help="Rebuild the index from scratch")
args = parser.parse_args()

# Ask question via input if not provided
if not args.query:
    try:
        question = input("\n❓ Please type your question: ").strip()
    except KeyboardInterrupt:
        print("\n❌ Aborted.")
        sys.exit(0)
else:
    question = args.query.strip()

if not question:
    print("⚠️ No question provided. Use default question.")
    #sys.exit(0)
    question = question_default

# Embedding and vector store setup
embed_model = HuggingFaceEmbedding(model_name="BAAI/bge-small-en-v1.5")
chroma_client = chromadb.PersistentClient(path=str(indexdir))
chroma_collection = chroma_client.get_or_create_collection("sfppy_kb")
vector_store = ChromaVectorStore(chroma_collection=chroma_collection)
storage_context = StorageContext.from_defaults(vector_store=vector_store)

# Build or load index
force_rebuild = True  # 🔁 toggle this flag manually
if force_rebuild or args.rebuild or not (indexdir / "chroma-collections.parquet").exists():
    print("📂 Building KB index from Markdown files...")
    documents = SimpleDirectoryReader(input_dir=str(kbdir), recursive=True).load_data()
    print(f"🔍 Indexed documents: {len(documents)}")
    index = VectorStoreIndex.from_documents(documents, storage_context=storage_context, embed_model=embed_model)
    storage_context.persist()
    print("🔄 Index persisted to:", indexdir)
else:
    print("✨ Loading prebuilt KB index...")
    index = load_index_from_storage(storage_context)

# Query engine setup
llm = Ollama(model="mistral", request_timeout=60.0, max_tokens=1024)
query_engine = index.as_query_engine(llm=llm, similarity_top_k=3)

# Run query
print("\n❓ Question:", question)
start = time.time()
response = query_engine.query(question)
end = time.time()

print("\n✅ Response:")
print(response.response)

print("\n⚠️ Context sources:")
for node in response.source_nodes:
    print(" -", node.metadata.get("file_path", "unknown"))

print(f"\n⏱️ Elapsed: {end - start:.2f} sec")



# %% Minimalist version (for those interested)
"""
    from llama_index.core import SimpleDirectoryReader, VectorStoreIndex
    from llama_index.embeddings.huggingface import HuggingFaceEmbedding
    from llama_index.llms.ollama import Ollama
    # Load KB
    documents = SimpleDirectoryReader(input_dir="./docs/KB", recursive=True).load_data()
    # Use a free HuggingFace embedder (no API key)
    embed_model = HuggingFaceEmbedding(model_name="BAAI/bge-small-en-v1.5")
    # Index documents
    index = VectorStoreIndex.from_documents(documents, embed_model=embed_model)
    # Define local LLM via Ollama
    llm = Ollama(model="mistral")
    # Query engine (RAG happens here)
    query_engine = index.as_query_engine(llm=llm)
    # Ask a question
    response = query_engine.query("What is the definition of a functional barrier in EU Regulation 10/2011?")
    print(response)
"""
