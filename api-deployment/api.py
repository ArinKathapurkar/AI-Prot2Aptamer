from flask import Flask, request, jsonify
from flask_cors import CORS
import pickle
import numpy as np
import os

app = Flask(__name__)
CORS(app)

# Load models at startup
model = pickle.load(open('my_aptamer_predictor.pkl', 'rb'))
kmer_library = pickle.load(open('kmer_library.pkl', 'rb'))
transition_matrix = pickle.load(open('transition_matrix.pkl', 'rb'))

@app.route('/api/generate', methods=['POST'])
def generate():
    data = request.json
    protein_seq = data['protein_sequence']
    
    # Your generation logic here
    # (simplified version - just the essentials)
    
    results = [
        {
            'rank': 1,
            'sequence': 'TCGTGGAGCTCC...',
            'length': 77,
            'gc_content': 66.2,
            'quality_score': 2.89
        }
        # ... generate real results
    ]
    
    return jsonify({'aptamers': results})

@app.route('/health', methods=['GET'])
def health():
    return jsonify({'status': 'healthy'})

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port)