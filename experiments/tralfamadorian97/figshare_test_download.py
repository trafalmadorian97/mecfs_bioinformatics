
import subprocess

def go():

    # Define your target file and authentication token
    file_url = "https://ndownloader.figshare.com/files/34517828"
    # file_url = "https://figshare.com/ndownloader/files/34517828"
    api_token = "b8f0ce9a33920dfd00857233cdb569578c25fec9401085e5489252bb26aed37e1991a3667fb7854ec3c1beee0448613e237532864f719942fa95ef4c650fbf39"

    # Build the aria2c command as a list of arguments
    aria2_cmd = [
        "pixi",
        "r",
        '--environment', 'download-env',
        "aria2c",
        # Pass the API token to bypass the WAF challenge
        # f"--header=Authorization: token {api_token}",
        # It is still good practice to spoof a standard User-Agent
        "--user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",
        # Set standard helpful aria2 flags for build systems (optional)
        "--quiet=false",
        "-l",
        "-",
        "--summary-interval=10",
        file_url
    ]
    print(aria2_cmd)

    print(f"Starting API-authenticated download for {file_url}...")

    # Execute the command
    try:
        subprocess.run(aria2_cmd, check=True)
        print("Download completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Download failed with exit code {e.returncode}")


if __name__ == "__main__":
    go()